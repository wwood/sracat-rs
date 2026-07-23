//! Command-line interface tests, exercising the built binary against a small
//! single-end fixture in tests/data.

use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use assert_cli::Assert;

fn fixture() -> String {
    format!(
        "{}/tests/data/ERR015558.lite.sra",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// Run the built binary with a wall-clock limit, so a hang (e.g. a writer/worker
/// deadlock) becomes a test failure instead of blocking the whole suite. Returns
/// whether the process exited successfully.
fn run_within(secs: u64, args: &[&str]) -> bool {
    let mut child = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn sracat-rs");
    let start = Instant::now();
    loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            return status.success();
        }
        if start.elapsed() > Duration::from_secs(secs) {
            let _ = child.kill();
            let _ = child.wait();
            panic!("sracat-rs blocked: no exit within {secs}s (args: {args:?})");
        }
        std::thread::sleep(Duration::from_millis(50));
    }
}

#[test]
fn help_succeeds() {
    Assert::main_binary()
        .with_args(&["--help"])
        .succeeds()
        .stdout()
        .contains("--single")
        .unwrap();
}

#[test]
fn missing_file_fails() {
    Assert::main_binary()
        .with_args(&["/no/such/file.sra"])
        .fails()
        .unwrap();
}

/// Aligned (cSRA) runs are extracted by default (READ reconstructed from the
/// alignment table). Uses the local fixture pulled by `pixi run fetch-testdata`
/// (ERR1540848: a small S. pneumoniae cSRA run with a PRIMARY_ALIGNMENT table).
/// Skips if the fixture has not been fetched.
#[test]
fn aligned_run_extracts_by_default() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping aligned_run_extracts_by_default: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    Assert::main_binary()
        .with_args(&["--single-out", "/dev/null", f.as_str()])
        .succeeds()
        .stdout()
        .contains(">ERR1540848.")
        .unwrap();
}

/// With --croak-on-aligned, the same cSRA run is refused rather than extracted.
/// Uses the fetched fixture; skips if absent.
#[test]
fn croak_on_aligned_refuses() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping croak_on_aligned_refuses: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    Assert::main_binary()
        .with_args(&["--croak-on-aligned", f.as_str()])
        .fails()
        .stderr()
        .contains("aligned")
        .unwrap();
}

/// Aligned reconstruction does not parallelise, so -t > 1 on an aligned run is
/// capped to a single thread (with a note) rather than running the pathological
/// multi-cursor path. Uses the fetched fixture; skips if absent.
#[test]
fn aligned_run_caps_threads_to_one() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping aligned_run_caps_threads_to_one: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    Assert::main_binary()
        .with_args(&["-t", "8", "--single-out", "/dev/null", f.as_str()])
        .succeeds()
        .stderr()
        .contains("single-threaded")
        .unwrap();
}

/// `-1`/`-2` split paired output: forward reads (with `/1`) to the first file,
/// reverse reads (with `/2`) to the second. Uses the aligned (paired) fixture;
/// skips if absent.
#[test]
fn split_output_to_read1_read2() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping split_output_to_read1_read2: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    let dir = std::env::temp_dir();
    let r1 = format!("{}/sracat_rs_split_r1.fasta", dir.display());
    let r2 = format!("{}/sracat_rs_split_r2.fasta", dir.display());
    let _ = std::fs::remove_file(&r1);
    let _ = std::fs::remove_file(&r2);
    Assert::main_binary()
        .with_args(&["-1", r1.as_str(), "-2", r2.as_str(), f.as_str()])
        .succeeds()
        .unwrap();
    let b1 = std::fs::read_to_string(&r1).expect("read1 written");
    let b2 = std::fs::read_to_string(&r2).expect("read2 written");
    assert!(
        b1.contains(">ERR1540848.") && b1.contains("/1"),
        "read1 has /1"
    );
    assert!(
        b2.contains(">ERR1540848.") && b2.contains("/2"),
        "read2 has /2"
    );
    assert!(!b1.contains("/2"), "read1 must not contain reverse reads");
    assert!(!b2.contains("/1"), "read2 must not contain forward reads");
    let _ = std::fs::remove_file(&r1);
    let _ = std::fs::remove_file(&r2);
}

/// Split pairs from a real unaligned paired run at both -t1 and -t4, and check
/// the parallel output is byte-identical to single-threaded. Uses DRR033172
/// (46,282 pairs over several decode chunks), fetched by `pixi run
/// fetch-testdata`; skips if absent. Timeout-guarded against a parallel hang.
#[test]
fn split_unaligned_pairs_t1_eq_t4() {
    let f = format!(
        "{}/tests/data/DRR033172/DRR033172.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping split_unaligned_pairs_t1_eq_t4: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    let dir = std::env::temp_dir();
    let p = |tag: &str| format!("{}/sracat_rs_{tag}.fasta", dir.display());
    let (a1, a2, b1, b2) = (p("u_t1_r1"), p("u_t1_r2"), p("u_t4_r1"), p("u_t4_r2"));
    for x in [&a1, &a2, &b1, &b2] {
        let _ = std::fs::remove_file(x);
    }

    assert!(
        run_within(60, &["-1", &a1, "-2", &a2, f.as_str()]),
        "t1 split did not succeed"
    );
    assert!(
        run_within(60, &["-t", "4", "-1", &b1, "-2", &b2, f.as_str()]),
        "t4 split did not succeed (blocked?)"
    );

    let read = |path: &str| std::fs::read_to_string(path).expect("split output written");
    let (r1a, r2a, r1b, r2b) = (read(&a1), read(&a2), read(&b1), read(&b2));

    // Parallel output must be byte-identical to single-threaded.
    assert_eq!(r1a, r1b, "read1 differs between -t1 and -t4");
    assert_eq!(r2a, r2b, "read2 differs between -t1 and -t4");

    // Correct split: equal counts, spanning >1 chunk, forward /1 and reverse /2.
    let heads = |s: &str| -> Vec<String> {
        s.lines()
            .filter(|l| l.starts_with('>'))
            .map(str::to_string)
            .collect()
    };
    let (h1, h2) = (heads(&r1a), heads(&r2a));
    assert_eq!(h1.len(), h2.len(), "unequal read1/read2 counts");
    assert!(
        h1.len() > 8192,
        "fixture should span multiple decode chunks"
    );
    assert!(
        h1.iter().all(|h| h.ends_with("/1")),
        "read1 headers must be /1"
    );
    assert!(
        h2.iter().all(|h| h.ends_with("/2")),
        "read2 headers must be /2"
    );

    for x in [&a1, &a2, &b1, &b2] {
        let _ = std::fs::remove_file(x);
    }
}

/// `-1`/`-2` with multiple threads must run the parallel writer to completion,
/// not deadlock. The single-end fixture exercises the parallel decode path; its
/// reads are orphans, so they go to --single-out. Timeout-guarded so a hang
/// fails the test rather than blocking CI forever.
#[test]
fn split_parallel_does_not_block() {
    let dir = std::env::temp_dir();
    let r1 = format!("{}/sracat_rs_par_r1.fasta", dir.display());
    let r2 = format!("{}/sracat_rs_par_r2.fasta", dir.display());
    let s = format!("{}/sracat_rs_par_s.fasta", dir.display());
    for f in [&r1, &r2, &s] {
        let _ = std::fs::remove_file(f);
    }
    let ok = run_within(
        30,
        &[
            "-t",
            "4",
            "-1",
            &r1,
            "-2",
            &r2,
            "--single-out",
            &s,
            fixture().as_str(),
        ],
    );
    assert!(ok, "parallel -1/-2 run did not exit successfully");
    let singles = std::fs::read_to_string(&s).expect("singles written");
    assert!(
        singles.contains("ERR015558"),
        "expected orphan reads routed to --single-out"
    );
    for f in [&r1, &r2, &s] {
        let _ = std::fs::remove_file(f);
    }
}

/// Split output to FIFOs must not block while opening the first mate FIFO before
/// the second. `--bench-read-only` avoids needing external readers; the test is
/// specifically about FIFO open behavior.
#[cfg(unix)]
#[test]
fn split_fifos_open_without_readers() {
    let dir = std::env::temp_dir();
    let r1 = format!("{}/sracat_rs_fifo_r1", dir.display());
    let r2 = format!("{}/sracat_rs_fifo_r2", dir.display());
    for f in [&r1, &r2] {
        let _ = std::fs::remove_file(f);
        let status = Command::new("mkfifo").arg(f).status().expect("run mkfifo");
        assert!(status.success(), "mkfifo failed for {f}");
    }

    let ok = run_within(
        30,
        &[
            "-1",
            &r1,
            "-2",
            &r2,
            "--bench-read-only",
            fixture().as_str(),
        ],
    );
    assert!(ok, "split FIFO outputs blocked while opening");

    for f in [&r1, &r2] {
        let _ = std::fs::remove_file(f);
    }
}

/// A real split-FIFO consumer can prefetch a large chunk from R1 before opening
/// R2. This mirrors `weebill`/needletail's initial R1 buffer fill and verifies
/// the CLI can keep writing until the delayed R2 reader catches up.
#[cfg(unix)]
#[test]
fn split_fifos_tolerate_delayed_r2_reader() {
    use std::io::Read;

    let f = format!(
        "{}/tests/data/DRR033172/DRR033172.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping split_fifos_tolerate_delayed_r2_reader: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }

    let dir = std::env::temp_dir();
    let tag = format!("{}_{}", std::process::id(), unique_nanos());
    let r1 = format!("{}/sracat_rs_delay_r1_{tag}", dir.display());
    let r2 = format!("{}/sracat_rs_delay_r2_{tag}", dir.display());
    for fifo in [&r1, &r2] {
        let status = Command::new("mkfifo")
            .arg(fifo)
            .status()
            .expect("run mkfifo");
        assert!(status.success(), "mkfifo failed for {fifo}");
    }

    let mut child = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(["-1", &r1, "-2", &r2, f.as_str()])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn sracat-rs");

    let (tx, rx) = std::sync::mpsc::channel();
    let r1_reader = {
        let r1 = r1.clone();
        std::thread::spawn(move || {
            let mut file = std::fs::File::open(&r1).expect("open r1 fifo");
            let mut prefetch = vec![0u8; 128 * 1024];
            file.read_exact(&mut prefetch).expect("prefetch r1");
            tx.send(()).expect("signal r1 prefetch");
            let mut rest = Vec::new();
            file.read_to_end(&mut rest).expect("drain r1");
            prefetch.len() + rest.len()
        })
    };

    rx.recv_timeout(Duration::from_secs(15))
        .expect("r1 prefetch should complete before r2 opens");

    let r2_reader = {
        let r2 = r2.clone();
        std::thread::spawn(move || {
            let mut file = std::fs::File::open(&r2).expect("open r2 fifo");
            let mut data = Vec::new();
            file.read_to_end(&mut data).expect("drain r2");
            data.len()
        })
    };

    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            break status;
        }
        if start.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            panic!("sracat-rs blocked with delayed r2 FIFO reader");
        }
        std::thread::sleep(Duration::from_millis(50));
    };
    assert!(status.success(), "sracat-rs did not exit successfully");

    let r1_bytes = r1_reader.join().expect("join r1 reader");
    let r2_bytes = r2_reader.join().expect("join r2 reader");
    assert!(r1_bytes > 128 * 1024, "expected substantial r1 output");
    assert!(r2_bytes > 128 * 1024, "expected substantial r2 output");

    for fifo in [&r1, &r2] {
        let _ = std::fs::remove_file(fifo);
    }
}

#[cfg(unix)]
fn unique_nanos() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos()
}

/// `-1` requires `-2` (and vice versa): clap rejects one without the other.
#[test]
fn read1_requires_read2() {
    Assert::main_binary()
        .with_args(&["-1", "/tmp/sracat_rs_only_r1.fasta", fixture().as_str()])
        .fails()
        .unwrap();
}

/// With `-1`/`-2` and an unpaired read but no single destination, refuse rather
/// than drop it (the fixture is single-end).
#[test]
fn split_without_single_sink_croaks() {
    Assert::main_binary()
        .with_args(&[
            "-1",
            "/tmp/sracat_rs_split_a.fasta",
            "-2",
            "/tmp/sracat_rs_split_b.fasta",
            fixture().as_str(),
        ])
        .fails()
        .stderr()
        .contains("unpaired")
        .unwrap();
}

#[test]
fn single_end_without_sink_croaks() {
    // The fixture is single-end, so with pairs going to stdout and no single
    // destination the tool must refuse rather than drop reads.
    Assert::main_binary()
        .with_args(&[fixture().as_str()])
        .fails()
        .stderr()
        .contains("unpaired")
        .unwrap();
}

#[test]
fn single_end_with_single_out_succeeds() {
    let out = format!(
        "{}/sracat_rs_test_single.fasta",
        std::env::temp_dir().display()
    );
    let _ = std::fs::remove_file(&out);
    Assert::main_binary()
        .with_args(&["--single-out", out.as_str(), fixture().as_str()])
        .succeeds()
        .unwrap();
    let body = std::fs::read_to_string(&out).expect("single output written");
    assert!(body.starts_with('>'), "expected FASTA output");
    assert!(
        body.contains("ERR015558"),
        "expected run name in read headers"
    );
    let _ = std::fs::remove_file(&out);
}

/// The single/orphan sink is opened lazily, so a cleanly-paired run with zero orphans
/// leaves no file. --eager-open-output opens it up front instead, so it exists (empty)
/// even with no orphans. Uses the paired DRR033172 fixture (46,282 pairs, 0 orphans);
/// skips if absent.
#[test]
fn eager_open_output_creates_empty_single_file() {
    let f = format!(
        "{}/tests/data/DRR033172/DRR033172.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping eager_open_output_creates_empty_single_file: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let lazy = format!("{}/sracat_rs_lazy_single_{pid}.fasta", dir.display());
    let eager = format!("{}/sracat_rs_eager_single_{pid}.fasta", dir.display());
    for p in [&lazy, &eager] {
        let _ = std::fs::remove_file(p);
    }

    // Default (lazy): a zero-orphan run never opens the sink, so no file is created.
    Assert::main_binary()
        .with_args(&["--single", lazy.as_str(), f.as_str()])
        .succeeds()
        .unwrap();
    assert!(
        !std::path::Path::new(&lazy).exists(),
        "lazy single sink must not be created for a zero-orphan run"
    );

    // --eager-open-output: the sink is opened at startup, so it exists but is empty.
    Assert::main_binary()
        .with_args(&[
            "--eager-open-output",
            "--single",
            eager.as_str(),
            f.as_str(),
        ])
        .succeeds()
        .unwrap();
    let meta = std::fs::metadata(&eager).expect("eager single sink must be created up front");
    assert_eq!(
        meta.len(),
        0,
        "a zero-orphan run must leave the eager single sink empty"
    );

    for p in [&lazy, &eager] {
        let _ = std::fs::remove_file(p);
    }
}

/// Regression for the streaming hang: when the single/orphan sink is a FIFO and the run
/// has zero orphans, the lazy open never happens so the pipe never gets a writer, and a
/// consumer blocking on open(O_RDONLY) hangs forever. --eager-open-output opens the pipe
/// up front, so the consumer instead sees a clean EOF. Timeout-guarded so a regression
/// fails the test rather than blocking forever.
#[cfg(unix)]
#[test]
fn eager_open_output_single_fifo_gets_clean_eof() {
    use std::io::Read;

    let f = format!(
        "{}/tests/data/DRR033172/DRR033172.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping eager_open_output_single_fifo_gets_clean_eof: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }

    let dir = std::env::temp_dir();
    let tag = format!("{}_{}", std::process::id(), unique_nanos());
    let s = format!("{}/sracat_rs_eager_fifo_{tag}", dir.display());
    let status = Command::new("mkfifo").arg(&s).status().expect("run mkfifo");
    assert!(status.success(), "mkfifo failed for {s}");

    // Pairs go to /dev/null; singles stream through the FIFO. Eager-open makes sracat
    // open the FIFO at startup (blocking until the reader below connects).
    let mut child = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(["--eager-open-output", "--single", &s, f.as_str()])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn sracat-rs");

    let s_reader = {
        let s = s.clone();
        std::thread::spawn(move || {
            let mut file = std::fs::File::open(&s).expect("open single fifo");
            let mut data = Vec::new();
            file.read_to_end(&mut data).expect("drain single fifo");
            data.len()
        })
    };

    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            break status;
        }
        if start.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            panic!("sracat-rs blocked streaming singles through a FIFO with --eager-open-output");
        }
        std::thread::sleep(Duration::from_millis(50));
    };
    assert!(status.success(), "sracat-rs did not exit successfully");

    let n = s_reader.join().expect("join single fifo reader");
    assert_eq!(n, 0, "a zero-orphan run must stream an empty singles pipe");

    let _ = std::fs::remove_file(&s);
}

/// --accept-singles streams single/orphan reads to stdout (interleaved with any
/// pairs) instead of demanding a separate destination. The fixture is single-end,
/// so with the default this croaks; with --accept-singles the reads land on
/// stdout.
#[test]
fn accept_singles_streams_singles_to_stdout() {
    Assert::main_binary()
        .with_args(&["--accept-singles", fixture().as_str()])
        .succeeds()
        .stdout()
        .contains(">ERR015558.")
        .unwrap();
}

/// --expect-singles mirrors the default: single/orphan reads go to stdout and a
/// paired spot would croak. The fixture is entirely single-end, so it succeeds
/// and the reads appear on stdout.
#[test]
fn expect_singles_streams_singles_to_stdout() {
    Assert::main_binary()
        .with_args(&["--expect-singles", fixture().as_str()])
        .succeeds()
        .stdout()
        .contains(">ERR015558.")
        .unwrap();
}

/// --sample N draws N spots at random (reading only those rows). For the
/// single-end fixture, N spots == N reads, so a sample smaller than the run
/// yields exactly N FASTA records, and the same --seed reproduces byte-identical
/// output. The sample size is calibrated to the fixture's read count (learned
/// from a full extraction) so the test holds whatever the fixture's size.
#[test]
fn sample_is_reproducible_and_bounded() {
    let run = |args: &[&str]| {
        let out = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
            .args(args)
            .output()
            .expect("run --sample");
        assert!(out.status.success(), "sample run failed: {args:?}");
        out.stdout
    };
    let count_reads = |b: &[u8]| {
        String::from_utf8_lossy(b)
            .lines()
            .filter(|l| l.starts_with('>'))
            .count()
    };

    // How many reads does the whole (single-end) run hold? Sample fewer than that
    // so we exercise the random-subset path, not the "extract everything" path.
    let total = count_reads(&run(&["--accept-singles", fixture().as_str()]));
    assert!(total >= 2, "fixture too small to subsample");
    let n = (total / 2).max(1).to_string();

    let a = run(&[
        "--accept-singles",
        "--sample",
        &n,
        "--seed",
        "42",
        fixture().as_str(),
    ]);
    let b = run(&[
        "--accept-singles",
        "--sample",
        &n,
        "--seed",
        "42",
        fixture().as_str(),
    ]);
    assert_eq!(a, b, "same seed must reproduce the same sample");

    let text = String::from_utf8(a).expect("utf8 output");
    let heads: Vec<&str> = text.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(
        heads.len(),
        n.parse::<usize>().unwrap(),
        "single-end: N sampled spots == N reads"
    );
    assert!(
        heads.iter().all(|h| h.contains("ERR015558")),
        "sampled headers carry the run name"
    );

    // The sample is emitted in random order, not sorted by row. Parse the row id
    // from each `>ERR015558.<row>` header and check they are not ascending. Only
    // assert when the sample is large enough that a shuffle landing in sorted
    // order by chance (probability 1/n!) is negligible.
    let rows: Vec<u64> = heads
        .iter()
        .map(|h| {
            h.rsplit('.')
                .next()
                .and_then(|r| {
                    r.trim_end_matches(|c: char| !c.is_ascii_digit())
                        .parse()
                        .ok()
                })
                .expect("row id in header")
        })
        .collect();
    if rows.len() >= 5 {
        assert!(
            rows.windows(2).any(|w| w[0] > w[1]),
            "sampled reads must be in random order, not sorted by row: {rows:?}"
        );
    }

    // A different seed should pick different rows. Only assert this when the pool
    // is large enough that two independent subsets colliding is negligible (on a
    // tiny fixture it could happen by chance).
    if total >= 8 {
        let c = run(&[
            "--accept-singles",
            "--sample",
            &n,
            "--seed",
            "7",
            fixture().as_str(),
        ]);
        assert_ne!(
            text.into_bytes(),
            c,
            "a different seed should draw a different sample"
        );
    }
}

/// --sample with N >= the run's spot count extracts the whole run, in storage
/// order — byte-identical to an unsampled extraction.
#[test]
fn sample_larger_than_run_extracts_everything() {
    let full = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(["--accept-singles", fixture().as_str()])
        .output()
        .expect("full run");
    assert!(full.status.success(), "full run failed");
    let sampled = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args([
            "--accept-singles",
            "--sample",
            "100000000",
            fixture().as_str(),
        ])
        .output()
        .expect("oversized sample run");
    assert!(sampled.status.success(), "oversized sample run failed");
    assert!(!full.stdout.is_empty(), "fixture should produce reads");
    assert_eq!(
        full.stdout, sampled.stdout,
        "sampling more spots than exist must equal a full extraction"
    );
}

/// --accept-singles and --expect-singles are mutually exclusive.
#[test]
fn accept_and_expect_singles_conflict() {
    Assert::main_binary()
        .with_args(&["--accept-singles", "--expect-singles", fixture().as_str()])
        .fails()
        .unwrap();
}

/// --accept-singles routes singles through the interleaved paired stream, so the
/// parallel writer path must run to completion (not deadlock) with -t > 1.
#[test]
fn accept_singles_parallel_does_not_block() {
    let ok = run_within(30, &["-t", "4", "--accept-singles", fixture().as_str()]);
    assert!(
        ok,
        "parallel --accept-singles run did not exit successfully"
    );
}

/// The interleaved singles+pairs stdout stream must stay byte-identical across
/// thread counts, like every other output mode.
#[test]
fn accept_singles_output_is_thread_stable() {
    let out1 = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(["--accept-singles", fixture().as_str()])
        .output()
        .expect("run --accept-singles at t1");
    assert!(out1.status.success(), "t1 run failed");
    let out4 = Command::new(env!("CARGO_BIN_EXE_sracat-rs"))
        .args(["-t", "4", "--accept-singles", fixture().as_str()])
        .output()
        .expect("run --accept-singles at t4");
    assert!(out4.status.success(), "t4 run failed");
    assert!(!out1.stdout.is_empty(), "expected reads on stdout");
    assert_eq!(
        out1.stdout, out4.stdout,
        "--accept-singles output must be byte-identical across thread counts"
    );
}
