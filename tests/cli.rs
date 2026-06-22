//! Command-line interface tests, exercising the built binary against a small
//! single-end fixture in tests/data.

use assert_cli::Assert;

fn fixture() -> String {
    format!(
        "{}/tests/data/ERR015558.lite.sra",
        env!("CARGO_MANIFEST_DIR")
    )
}

#[test]
fn help_succeeds() {
    Assert::main_binary()
        .with_args(&["--help"])
        .succeeds()
        .stdout()
        .contains("--single-out")
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
