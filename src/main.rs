//! sracat-rs: fast, deterministic extraction of reads from SRA files.
//!
//! Reads the SEQUENCE table of an SRA run in storage (row) order via the
//! ncbi-vdb cursor API (through a small C shim). Output is therefore repeatable
//! across runs. Paired spots are emitted interleaved; single/orphan spots are
//! routed to a separate stream. Aligned runs (where READ is reconstructed from
//! alignments) are refused unless --allow-aligned is given.

use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
#[cfg(unix)]
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};

mod ffi;
use ffi::Run;

/// Extract reads from SRA files as FASTA/FASTQ, in repeatable storage order.
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    /// SRA files (.sra) or accessions.
    #[arg(required = true, value_name = "SRA")]
    inputs: Vec<String>,

    /// Write paired reads to <PREFIX>.paired.{fasta,fastq} and singles to
    /// <PREFIX>.single.{fasta,fastq} instead of streaming pairs to stdout.
    #[arg(short = 'o', long, value_name = "PREFIX")]
    output_prefix: Option<String>,

    /// Split paired output: write the forward read of each pair to this file
    /// (requires -2). Mates are not interleaved. Single/orphan reads still go to
    /// --single (or -o), or the run croaks if neither is given.
    #[arg(
        short = '1',
        long = "read1",
        value_name = "FILE",
        requires = "read2",
        conflicts_with = "output_prefix"
    )]
    read1: Option<String>,

    /// Split paired output: write the reverse read of each pair to this file
    /// (requires -1).
    #[arg(
        short = '2',
        long = "read2",
        value_name = "FILE",
        requires = "read1",
        conflicts_with = "output_prefix"
    )]
    read2: Option<String>,

    /// Write single/orphan reads to this file (when pairs are streamed to stdout
    /// or split via -1/-2).
    #[arg(long = "single", alias = "single-out", value_name = "FILE")]
    single_out: Option<String>,

    /// Open the single/orphan output (--single, or -o's *.single) up front instead
    /// of lazily on the first orphan read. By default the sink is opened only when
    /// an orphan actually appears, so a cleanly-paired run leaves no empty file --
    /// but when the destination is a FIFO/named pipe, that lazy open means the pipe
    /// never gets a writer for a run with no orphans, and a reader blocking on
    /// open(O_RDONLY) hangs forever waiting for a partner. Pass this when streaming
    /// singles through a FIFO so the reader always gets a clean EOF.
    #[arg(long)]
    eager_open_output: bool,

    /// Stream single/orphan reads to stdout interleaved with the paired reads,
    /// in storage order, instead of requiring a separate destination for them.
    /// All reads (pairs and singles) are written through the one stdout stream,
    /// so records stay intact.
    #[arg(
        long,
        conflicts_with_all = ["output_prefix", "read1", "read2", "single_out", "expect_singles"]
    )]
    accept_singles: bool,

    /// Invert the default single-handling: stream single/orphan reads to stdout
    /// and croak if any paired spot is encountered. This mirrors the default
    /// (which streams pairs to stdout and croaks on a single/orphan read).
    #[arg(
        long,
        conflicts_with_all = ["output_prefix", "read1", "read2", "single_out", "accept_singles"]
    )]
    expect_singles: bool,

    /// Write FASTQ (with quality scores) instead of FASTA.
    #[arg(long)]
    qual: bool,

    /// Include technical reads (default: biological reads only).
    #[arg(long)]
    include_technical: bool,

    /// Refuse aligned (cSRA) runs instead of extracting them. By default aligned
    /// runs are extracted by reading the computed READ column, which ncbi-vdb
    /// reconstructs per spot from the alignment table (correct, spot-ordered, no
    /// temp files) at the cost of random access into that table; pass this flag
    /// to error out on such runs instead.
    #[arg(long)]
    croak_on_aligned: bool,

    /// Number of extraction threads. >1 decodes contiguous row ranges in
    /// parallel and writes them in order through a single writer (still
    /// repeatable; no temp files).
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Show a per-input progress bar (spots processed) on stderr.
    #[arg(long)]
    progress: bool,

    /// (benchmark) read + classify spots but skip all formatting/output.
    #[arg(long, hide = true)]
    bench_read_only: bool,
}

#[derive(Clone, Copy)]
struct Opts {
    qual: bool,
    include_technical: bool,
    allow_aligned: bool,
    bench_read_only: bool,
    /// --accept-singles: route single/orphan reads into the interleaved paired
    /// stream (stdout) rather than the separate single destination.
    singles_to_paired: bool,
    /// --expect-singles: refuse any paired spot (only singles are wanted).
    croak_on_paired: bool,
}

#[derive(Default, Clone, Copy)]
struct Counts {
    pairs: u64,
    singles: u64,
    skipped: u64,
}

impl Counts {
    fn add(&mut self, o: Counts) {
        self.pairs += o.pairs;
        self.singles += o.singles;
        self.skipped += o.skipped;
    }
}

/// Where paired reads go. Either one interleaved stream (stdout or
/// `<prefix>.paired`), or two separate streams with the forward read in the
/// first and the reverse read in the second (`-1`/`-2`).
enum PairedSink<'a> {
    Interleaved(&'a mut dyn Write),
    Split(&'a mut dyn Write, &'a mut dyn Write),
}

impl PairedSink<'_> {
    /// True if mates are split across two streams (`-1`/`-2`).
    fn is_split(&self) -> bool {
        matches!(self, PairedSink::Split(..))
    }
}

/// Write-buffer size for FIFO outputs. Kept well under the enlarged 1 MiB pipe
/// capacity so a single writer keeps split R1/R2 pipes finely interleaved (see
/// `open_split_output`); large enough that per-field record writes coalesce into
/// a handful of `write()` syscalls per record instead of one syscall per field.
const FIFO_WRITE_BUF: usize = 64 << 10;

#[cfg(unix)]
fn open_split_output(path: &str) -> io::Result<Box<dyn Write>> {
    use std::os::unix::fs::FileTypeExt;

    match fs::metadata(path) {
        Ok(meta) if meta.file_type().is_fifo() => {
            // FIFO consumers may open/read split mates sequentially. For example,
            // needletail fills a 64 KiB buffer from R1 before `weebill` advances
            // to R2, so a default 64 KiB R2 pipe can deadlock even when sracat
            // alternates R1/R2 writes. O_RDWR lets both FIFO opens complete; a
            // larger pipe gives the delayed mate enough room until its reader
            // catches up.
            let f = OpenOptions::new().read(true).write(true).open(path)?;
            enlarge_pipe_buffer(&f);
            // Buffer FIFO writes: the single-threaded extractor emits each record
            // field-by-field, so an unbuffered FIFO turns every field into a tiny
            // write() syscall (~16 bytes each). That pins a core spinning in the
            // kernel and starves the reader into per-record blocking reads,
            // collapsing pipeline throughput (e.g. a 9 s sketch stretches past
            // 100 s). The buffer must stay well below the pipe capacity (enlarged
            // to 1 MiB above): with a split R1/R2 run, a single writer thread
            // flushes the two pipes in turn, so a flush as large as the pipe can
            // fill R1 completely while R2 is still buffered, deadlocking a reader
            // that needs R2 to advance. A small buffer keeps R1/R2 finely
            // interleaved so neither pipe fills while the other sits unflushed.
            Ok(Box::new(BufWriter::with_capacity(FIFO_WRITE_BUF, f)) as Box<dyn Write>)
        }
        _ => File::create(path)
            .map(|f| Box::new(BufWriter::with_capacity(1 << 20, f)) as Box<dyn Write>),
    }
}

#[cfg(target_os = "linux")]
fn enlarge_pipe_buffer(f: &File) {
    const FIFO_PIPE_BYTES: libc::c_int = 1 << 20;
    unsafe {
        // Best effort: unprivileged callers can usually raise small pipes up to
        // /proc/sys/fs/pipe-max-size. If the kernel refuses, ordinary FIFO
        // semantics still apply and the write will surface any real I/O error.
        libc::fcntl(f.as_raw_fd(), libc::F_SETPIPE_SZ, FIFO_PIPE_BYTES);
    }
}

#[cfg(all(unix, not(target_os = "linux")))]
fn enlarge_pipe_buffer(_: &File) {}

#[cfg(not(unix))]
fn open_split_output(path: &str) -> io::Result<Box<dyn Write>> {
    File::create(path).map(|f| Box::new(BufWriter::with_capacity(1 << 20, f)) as Box<dyn Write>)
}

/// Destination for single/orphan reads, opened lazily so no empty file is
/// created for a cleanly paired run, and so the "no destination" case can fail
/// only if a single read actually appears.
struct SingleWriter {
    dest: SingleDest,
    inner: Option<BufWriter<Box<dyn Write>>>,
}

enum SingleDest {
    /// No destination configured: refuse if any single read is written.
    Fail,
    Path(PathBuf),
    /// --expect-singles: single/orphan reads are the wanted output, streamed
    /// to stdout.
    Stdout,
}

impl SingleWriter {
    fn ensure(&mut self) -> io::Result<&mut BufWriter<Box<dyn Write>>> {
        if self.inner.is_none() {
            let w: Box<dyn Write> = match &self.dest {
                SingleDest::Fail => {
                    return Err(io::Error::other(
                        "encountered unpaired read(s) but no destination for them; \
                         pass --single <file>, -o <prefix>, or --accept-singles",
                    ))
                }
                SingleDest::Path(p) => Box::new(
                    File::create(p)
                        .map_err(|e| io::Error::other(format!("creating {}: {e}", p.display())))?,
                ),
                SingleDest::Stdout => Box::new(io::stdout().lock()),
            };
            self.inner = Some(BufWriter::with_capacity(1 << 20, w));
        }
        Ok(self.inner.as_mut().unwrap())
    }

    fn finish(mut self) -> io::Result<()> {
        if let Some(w) = self.inner.as_mut() {
            w.flush()?;
        }
        Ok(())
    }
}

impl Write for SingleWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.ensure()?.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        match self.inner.as_mut() {
            Some(w) => w.flush(),
            None => Ok(()),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let opts = Opts {
        qual: cli.qual,
        include_technical: cli.include_technical,
        // Aligned runs are extracted by default; --croak-on-aligned opts back
        // into refusing them.
        allow_aligned: !cli.croak_on_aligned,
        bench_read_only: cli.bench_read_only,
        singles_to_paired: cli.accept_singles,
        croak_on_paired: cli.expect_singles,
    };
    let ext = if cli.qual { "fastq" } else { "fasta" };

    // Paired output: split to two files (-1/-2), else interleaved to one stream
    // (<prefix>.paired or stdout). The owned writers live here; `paired` borrows
    // them and is dropped before they are flushed below.
    let mut inter: Option<BufWriter<Box<dyn Write>>> = None;
    let mut split: Option<(Box<dyn Write>, Box<dyn Write>)> = None;
    match (&cli.read1, &cli.read2) {
        (Some(r1), Some(r2)) => {
            let f1 = open_split_output(r1).with_context(|| format!("creating {r1}"))?;
            let f2 = open_split_output(r2).with_context(|| format!("creating {r2}"))?;
            split = Some((f1, f2));
        }
        // clap guarantees -1 and -2 come together, so the remaining cases have
        // neither set: fall back to an interleaved stream.
        _ => {
            inter = Some(match &cli.output_prefix {
                Some(prefix) => {
                    let path = format!("{prefix}.paired.{ext}");
                    let f = File::create(&path).with_context(|| format!("creating {path}"))?;
                    BufWriter::with_capacity(1 << 20, Box::new(f) as Box<dyn Write>)
                }
                // --expect-singles streams singles to stdout and refuses pairs, so
                // the paired sink is never written; discard into a sink rather than
                // hold a second handle on stdout alongside the single writer.
                None if cli.expect_singles => {
                    BufWriter::with_capacity(1 << 20, Box::new(io::sink()) as Box<dyn Write>)
                }
                None => BufWriter::with_capacity(1 << 20, Box::new(io::stdout().lock())),
            });
        }
    }
    let mut single = SingleWriter {
        dest: match (&cli.output_prefix, &cli.single_out, cli.expect_singles) {
            // --expect-singles: single/orphan reads are the wanted output.
            (_, _, true) => SingleDest::Stdout,
            (Some(prefix), _, _) => SingleDest::Path(format!("{prefix}.single.{ext}").into()),
            (None, Some(path), _) => SingleDest::Path(path.into()),
            // --accept-singles routes singles into the paired stream, so its
            // Fail destination is never reached.
            (None, None, _) => SingleDest::Fail,
        },
        inner: None,
    };

    // --eager-open-output: open the single sink now rather than on the first orphan.
    // For a FIFO destination this guarantees the pipe gets a writer even when the run
    // has no orphans, so a consumer reading the pipe sees a clean EOF instead of
    // blocking forever on open(O_RDONLY). A Fail/Stdout destination has nothing to
    // pre-open (Fail must stay lazy so it only errors on a real orphan).
    if cli.eager_open_output {
        if let SingleDest::Path(_) = single.dest {
            single
                .ensure()
                .context("eagerly opening the single/orphan output")?;
        }
    }

    let threads = cli.threads.max(1);
    let mut totals = Counts::default();
    // Scoped so `paired`'s borrow of inter/split ends before they are flushed.
    {
        let mut paired = match (&mut inter, &mut split) {
            (Some(w), _) => PairedSink::Interleaved(w),
            (_, Some((w1, w2))) => PairedSink::Split(w1, w2),
            _ => unreachable!("one of inter/split is always set"),
        };
        for input in &cli.inputs {
            let name = derive_name(input);
            // Open once up front. Aligned (cSRA) reconstruction does random
            // access into the alignment table that does not parallelise — with
            // multiple cursors it degrades catastrophically — so hard-cap
            // aligned runs to a single thread regardless of -t.
            let run = Run::open(input, opts.qual, opts.allow_aligned)?;
            let eff_threads = if run.is_aligned() { 1 } else { threads };
            if eff_threads < threads {
                eprintln!(
                    "note: {name} is aligned (cSRA); extracting single-threaded (-t{threads} ignored)"
                );
            }
            // Per-input progress bar (created after the note above so it doesn't
            // clobber that line). Workers increment it as they classify spots.
            let pb = cli.progress.then(|| make_progress(&name, run.row_count()));
            let c = if eff_threads == 1 {
                let (lo, hi) = (run.first_row(), run.first_row() + run.row_count() as i64);
                extract_range(
                    &run,
                    lo,
                    hi,
                    &name,
                    &mut paired,
                    &mut single,
                    opts,
                    pb.as_ref(),
                )?
            } else {
                drop(run); // workers open their own cursors in extract_parallel
                extract_parallel(
                    input,
                    &name,
                    eff_threads,
                    &mut paired,
                    &mut single,
                    opts,
                    pb.as_ref(),
                )?
            };
            if let Some(pb) = pb {
                pb.finish();
            }
            totals.add(c);
        }
    }

    if let Some(mut w) = inter {
        w.flush()?;
    }
    if let Some((mut w1, mut w2)) = split {
        w1.flush()?;
        w2.flush()?;
    }
    single.finish()?;
    eprintln!("spots paired   : {}", totals.pairs);
    eprintln!("reads single   : {}", totals.singles);
    if totals.skipped > 0 {
        eprintln!("spots skipped  : {} (no biological reads)", totals.skipped);
    }
    Ok(())
}

/// Build a per-input progress bar over `total` spots, labelled with the run
/// name and rendered on stderr.
fn make_progress(name: &str, total: u64) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::with_template(
            "{msg}: [{bar:40.cyan/blue}] {human_pos}/{human_len} spots ({percent}%) {per_sec} ETA {eta}",
        )
        .expect("static progress template")
        .progress_chars("=>-"),
    );
    pb.set_message(name.to_string());
    pb
}

/// Extract rows `[lo, hi)` of an opened run, writing paired reads to `paired`
/// (interleaved or split) and single reads to `single`. When `progress` is set,
/// the bar is advanced once per spot; increments are batched to keep the atomic
/// traffic low and are safe to call concurrently from parallel workers.
#[allow(clippy::too_many_arguments)]
fn extract_range(
    run: &Run,
    lo: i64,
    hi: i64,
    name: &str,
    paired: &mut PairedSink<'_>,
    single: &mut dyn Write,
    opts: Opts,
    progress: Option<&ProgressBar>,
) -> Result<Counts> {
    let mut counts = Counts::default();
    let mut sel: Vec<(usize, usize)> = Vec::new();
    let mut qbuf: Vec<u8> = Vec::new();
    // Batch progress ticks so a huge run isn't a storm of atomic increments.
    const TICK: u64 = 4096;
    let mut since_tick = 0u64;

    for row in lo..hi {
        let spot = run.read_spot(row)?;

        sel.clear();
        let mut off = 0usize;
        for (&len32, &ty) in spot.read_len.iter().zip(spot.read_type.iter()) {
            let len = len32 as usize;
            if opts.include_technical || (ty & 1) != 0 {
                sel.push((off, len));
            }
            off += len;
        }

        match sel.len() {
            0 => counts.skipped += 1,
            1 => {
                if !opts.bench_read_only {
                    let (o, l) = sel[0];
                    // --accept-singles: interleave the single read into the paired
                    // stdout stream (in row order) instead of the separate single
                    // destination, so pairs and singles share one intact stream.
                    if opts.singles_to_paired {
                        match paired {
                            PairedSink::Interleaved(w) => {
                                write_read(&mut **w, name, row, None, &spot, o, l, &mut qbuf)?
                            }
                            PairedSink::Split(..) => unreachable!(
                                "--accept-singles is incompatible with -1/-2 split output"
                            ),
                        }
                    } else {
                        write_read(single, name, row, None, &spot, o, l, &mut qbuf)?;
                    }
                }
                counts.singles += 1;
            }
            2 => {
                if opts.croak_on_paired {
                    bail!(
                        "{name}: spot {row} is a paired spot, but --expect-singles \
                         was given (only single/orphan reads are permitted)"
                    );
                }
                if !opts.bench_read_only {
                    let (o0, l0) = sel[0];
                    let (o1, l1) = sel[1];
                    match paired {
                        PairedSink::Interleaved(w) => {
                            write_read(&mut **w, name, row, Some(1), &spot, o0, l0, &mut qbuf)?;
                            write_read(&mut **w, name, row, Some(2), &spot, o1, l1, &mut qbuf)?;
                        }
                        PairedSink::Split(w1, w2) => {
                            write_read(&mut **w1, name, row, Some(1), &spot, o0, l0, &mut qbuf)?;
                            write_read(&mut **w2, name, row, Some(2), &spot, o1, l1, &mut qbuf)?;
                        }
                    }
                }
                counts.pairs += 1;
            }
            n => bail!(
                "{name}: spot {row} has {n} biological reads (>2); not supported \
                 (use --include-technical to inspect, or file an issue)"
            ),
        }
        if let Some(pb) = progress {
            since_tick += 1;
            if since_tick >= TICK {
                pb.inc(since_tick);
                since_tick = 0;
            }
        }
    }
    if let Some(pb) = progress {
        pb.inc(since_tick);
    }
    Ok(counts)
}

/// Decode the run in parallel and write in order. Worker threads each own a
/// cursor, pull contiguous row chunks via an atomic counter, format them into
/// in-memory buffers, and hand them to the writer (this thread), which emits
/// chunks in index order. Output is byte-identical to the single-threaded path.
/// Memory is bounded by a window on how far decoding may run ahead of writing.
fn extract_parallel(
    input: &str,
    name: &str,
    threads: usize,
    paired: &mut PairedSink<'_>,
    single: &mut dyn Write,
    opts: Opts,
    progress: Option<&ProgressBar>,
) -> Result<Counts> {
    use std::collections::BTreeMap;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::mpsc::sync_channel;

    // Validate (incl. aligned check) and learn the row range up front.
    let (first, count) = {
        let run = Run::open(input, opts.qual, opts.allow_aligned)?;
        (run.first_row(), run.row_count())
    };
    let hi = first + count as i64;
    let split = paired.is_split(); // whether workers format two buffers or one

    const CHUNK: i64 = 8192; // spots per work unit
    let window = threads as u64 * 4; // max chunks decoded ahead of the writer

    let next_chunk = AtomicU64::new(0); // next chunk index to hand out
    let next_write = AtomicU64::new(0); // next chunk index still to be written

    // Formatted output for one chunk: (forward/interleaved, reverse, single).
    type Bufs = (Vec<u8>, Vec<u8>, Vec<u8>);
    type Msg = Result<(u64, Bufs, Counts)>;
    let (tx, rx) = sync_channel::<Msg>(threads * 2);

    std::thread::scope(|scope| -> Result<Counts> {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next_chunk, next_write) = (&next_chunk, &next_write);
            scope.spawn(move || {
                let run = match Run::open(input, opts.qual, opts.allow_aligned) {
                    Ok(r) => r,
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                };
                loop {
                    let idx = next_chunk.load(Ordering::Acquire);
                    if first + idx as i64 * CHUNK >= hi {
                        break;
                    }
                    // Keep memory bounded: don't decode too far ahead of writing.
                    if idx >= next_write.load(Ordering::Acquire) + window {
                        std::thread::yield_now();
                        continue;
                    }
                    if next_chunk
                        .compare_exchange(idx, idx + 1, Ordering::AcqRel, Ordering::Acquire)
                        .is_err()
                    {
                        continue;
                    }
                    let lo = first + idx as i64 * CHUNK;
                    let chi = (lo + CHUNK).min(hi);
                    let mut p1 = Vec::new(); // forward, or interleaved when not split
                    let mut p2 = Vec::new(); // reverse (split only)
                    let mut sbuf = Vec::new();
                    let res = {
                        let mut ps = if split {
                            PairedSink::Split(&mut p1, &mut p2)
                        } else {
                            PairedSink::Interleaved(&mut p1)
                        };
                        extract_range(&run, lo, chi, name, &mut ps, &mut sbuf, opts, progress)
                    };
                    let msg = res.map(|c| (idx, (p1, p2, sbuf), c));
                    let failed = msg.is_err();
                    if tx.send(msg).is_err() || failed {
                        break;
                    }
                }
            });
        }
        drop(tx); // workers hold the only remaining senders

        // Writer: reorder by chunk index and emit consecutively.
        let mut counts = Counts::default();
        let mut pending: BTreeMap<u64, Bufs> = BTreeMap::new();
        let mut expected = 0u64;
        for msg in rx {
            let (idx, bufs, c) = msg?;
            counts.add(c);
            pending.insert(idx, bufs);
            while let Some((p1, p2, sbuf)) = pending.remove(&expected) {
                match &mut *paired {
                    PairedSink::Interleaved(w) => w.write_all(&p1)?,
                    PairedSink::Split(w1, w2) => write_split_interleaved(*w1, *w2, &p1, &p2)?,
                }
                if !sbuf.is_empty() {
                    single.write_all(&sbuf)?;
                }
                expected += 1;
                next_write.store(expected, Ordering::Release);
            }
        }
        Ok(counts)
    })
}

/// Write a decoded chunk's forward (`p1`) and reverse (`p2`) buffers to the two
/// split (`-1`/`-2`) outputs, interleaved in small blocks. The parallel writer
/// is a single thread feeding both streams; when the outputs are FIFOs, writing
/// a whole multi-megabyte chunk to `w1` before `w2` fills (and then blocks on)
/// the R1 pipe while the R2 pipe stays empty, deadlocking a reader that consumes
/// mates in lockstep. `write_all` of a slice larger than the FIFO's `BufWriter`
/// bypasses that buffer, so the buffer alone does not prevent it. Interleaving
/// bounds the R1/R2 skew to `BLK`, which stays well under the enlarged pipe
/// capacity. For regular files this is simply two buffered writers filling in
/// step, at no meaningful cost.
fn write_split_interleaved(
    w1: &mut dyn Write,
    w2: &mut dyn Write,
    p1: &[u8],
    p2: &[u8],
) -> io::Result<()> {
    const BLK: usize = 32 << 10;
    let (mut a, mut b) = (p1, p2);
    while !a.is_empty() || !b.is_empty() {
        if !a.is_empty() {
            let n = a.len().min(BLK);
            w1.write_all(&a[..n])?;
            a = &a[n..];
        }
        if !b.is_empty() {
            let n = b.len().min(BLK);
            w2.write_all(&b[..n])?;
            b = &b[n..];
        }
    }
    Ok(())
}

/// Write one read as FASTA or FASTQ. `(off, len)` selects the read within the
/// spot; `qbuf` is a reused scratch buffer for phred->ASCII conversion.
#[allow(clippy::too_many_arguments)]
fn write_read(
    w: &mut dyn Write,
    name: &str,
    row: i64,
    mate: Option<u8>,
    spot: &ffi::Spot<'_>,
    off: usize,
    len: usize,
    qbuf: &mut Vec<u8>,
) -> io::Result<()> {
    let seq = &spot.bases[off..off + len];
    match spot.quals {
        Some(all_q) => {
            match mate {
                Some(m) => writeln!(w, "@{name}.{row}/{m}")?,
                None => writeln!(w, "@{name}.{row}")?,
            }
            w.write_all(seq)?;
            w.write_all(b"\n+\n")?;
            qbuf.clear();
            qbuf.extend(all_q[off..off + len].iter().map(|&q| q.saturating_add(33)));
            w.write_all(qbuf)?;
            w.write_all(b"\n")
        }
        None => {
            match mate {
                Some(m) => writeln!(w, ">{name}.{row}/{m}")?,
                None => writeln!(w, ">{name}.{row}")?,
            }
            w.write_all(seq)?;
            w.write_all(b"\n")
        }
    }
}

/// Run name from an input path: basename with a trailing `.sra` stripped.
fn derive_name(input: &str) -> String {
    let p = Path::new(input);
    let is_sra = p.extension().is_some_and(|e| e.eq_ignore_ascii_case("sra"));
    let stem = if is_sra { p.file_stem() } else { p.file_name() };
    stem.and_then(|s| s.to_str()).unwrap_or(input).to_string()
}
