//! sracat-rs: fast, deterministic extraction of reads from SRA files.
//!
//! Reads the SEQUENCE table of an SRA run in storage (row) order via the
//! ncbi-vdb cursor API (through a small C shim). Output is therefore repeatable
//! across runs. Paired spots are emitted interleaved; single/orphan spots are
//! routed to a separate stream. Aligned runs (where READ is reconstructed from
//! alignments) are refused.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Parser;

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

    /// When streaming pairs to stdout, write single/orphan reads to this file.
    #[arg(long, value_name = "FILE")]
    single_out: Option<String>,

    /// Write FASTQ (with quality scores) instead of FASTA.
    #[arg(long)]
    qual: bool,

    /// Include technical reads (default: biological reads only).
    #[arg(long)]
    include_technical: bool,

    /// Number of extraction threads. >1 decodes contiguous row ranges in
    /// parallel and writes them in order through a single writer (still
    /// repeatable; no temp files).
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// (benchmark) read + classify spots but skip all formatting/output.
    #[arg(long, hide = true)]
    bench_read_only: bool,
}

#[derive(Clone, Copy)]
struct Opts {
    qual: bool,
    include_technical: bool,
    bench_read_only: bool,
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

/// Destination for single/orphan reads, opened lazily so no empty file is
/// created for a cleanly paired run, and so the "no destination" case can fail
/// only if a single read actually appears.
struct SingleWriter {
    dest: SingleDest,
    inner: Option<BufWriter<File>>,
}

enum SingleDest {
    /// No destination configured: refuse if any single read is written.
    Fail,
    Path(PathBuf),
}

impl SingleWriter {
    fn ensure(&mut self) -> io::Result<&mut BufWriter<File>> {
        if self.inner.is_none() {
            match &self.dest {
                SingleDest::Fail => {
                    return Err(io::Error::other(
                        "encountered unpaired read(s) but no destination for them; \
                         pass --single-out <file> or -o <prefix>",
                    ))
                }
                SingleDest::Path(p) => {
                    let f = File::create(p)
                        .map_err(|e| io::Error::other(format!("creating {}: {e}", p.display())))?;
                    self.inner = Some(BufWriter::with_capacity(1 << 20, f));
                }
            }
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
        bench_read_only: cli.bench_read_only,
    };
    let ext = if cli.qual { "fastq" } else { "fasta" };

    let mut paired: BufWriter<Box<dyn Write>> = match &cli.output_prefix {
        Some(prefix) => {
            let path = format!("{prefix}.paired.{ext}");
            let f = File::create(&path).with_context(|| format!("creating {path}"))?;
            BufWriter::with_capacity(1 << 20, Box::new(f))
        }
        None => BufWriter::with_capacity(1 << 20, Box::new(io::stdout().lock())),
    };

    let mut single = SingleWriter {
        dest: match (&cli.output_prefix, &cli.single_out) {
            (Some(prefix), _) => SingleDest::Path(format!("{prefix}.single.{ext}").into()),
            (None, Some(path)) => SingleDest::Path(path.into()),
            (None, None) => SingleDest::Fail,
        },
        inner: None,
    };

    let threads = cli.threads.max(1);
    let mut totals = Counts::default();
    for input in &cli.inputs {
        let name = derive_name(input);
        let c = if threads == 1 {
            let run = Run::open(input, opts.qual)?;
            let (lo, hi) = (run.first_row(), run.first_row() + run.row_count() as i64);
            extract_range(&run, lo, hi, &name, &mut paired, &mut single, opts)?
        } else {
            extract_parallel(input, &name, threads, &mut paired, &mut single, opts)?
        };
        totals.add(c);
    }

    paired.flush()?;
    single.finish()?;
    eprintln!("spots paired   : {}", totals.pairs);
    eprintln!("reads single   : {}", totals.singles);
    if totals.skipped > 0 {
        eprintln!("spots skipped  : {} (no biological reads)", totals.skipped);
    }
    Ok(())
}

/// Extract rows `[lo, hi)` of an opened run, writing paired reads to `paired`
/// and single reads to `single`.
fn extract_range(
    run: &Run,
    lo: i64,
    hi: i64,
    name: &str,
    paired: &mut dyn Write,
    single: &mut dyn Write,
    opts: Opts,
) -> Result<Counts> {
    let mut counts = Counts::default();
    let mut sel: Vec<(usize, usize)> = Vec::new();
    let mut qbuf: Vec<u8> = Vec::new();

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
                    write_read(single, name, row, None, &spot, o, l, &mut qbuf)?;
                }
                counts.singles += 1;
            }
            2 => {
                if !opts.bench_read_only {
                    let (o0, l0) = sel[0];
                    let (o1, l1) = sel[1];
                    write_read(paired, name, row, Some(1), &spot, o0, l0, &mut qbuf)?;
                    write_read(paired, name, row, Some(2), &spot, o1, l1, &mut qbuf)?;
                }
                counts.pairs += 1;
            }
            n => bail!(
                "{name}: spot {row} has {n} biological reads (>2); not supported \
                 (use --include-technical to inspect, or file an issue)"
            ),
        }
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
    paired: &mut dyn Write,
    single: &mut dyn Write,
    opts: Opts,
) -> Result<Counts> {
    use std::collections::BTreeMap;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::mpsc::sync_channel;

    // Validate (incl. aligned check) and learn the row range up front.
    let (first, count) = {
        let run = Run::open(input, opts.qual)?;
        (run.first_row(), run.row_count())
    };
    let hi = first + count as i64;

    const CHUNK: i64 = 8192; // spots per work unit
    let window = threads as u64 * 4; // max chunks decoded ahead of the writer

    let next_chunk = AtomicU64::new(0); // next chunk index to hand out
    let next_write = AtomicU64::new(0); // next chunk index still to be written

    type Msg = Result<(u64, Vec<u8>, Vec<u8>, Counts)>;
    let (tx, rx) = sync_channel::<Msg>(threads * 2);

    std::thread::scope(|scope| -> Result<Counts> {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next_chunk, next_write) = (&next_chunk, &next_write);
            scope.spawn(move || {
                let run = match Run::open(input, opts.qual) {
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
                    let mut pbuf = Vec::new();
                    let mut sbuf = Vec::new();
                    let msg = extract_range(&run, lo, chi, name, &mut pbuf, &mut sbuf, opts)
                        .map(|c| (idx, pbuf, sbuf, c));
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
        let mut pending: BTreeMap<u64, (Vec<u8>, Vec<u8>)> = BTreeMap::new();
        let mut expected = 0u64;
        for msg in rx {
            let (idx, pbuf, sbuf, c) = msg?;
            counts.add(c);
            pending.insert(idx, (pbuf, sbuf));
            while let Some((pbuf, sbuf)) = pending.remove(&expected) {
                paired.write_all(&pbuf)?;
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
