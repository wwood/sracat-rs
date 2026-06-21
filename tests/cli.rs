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

/// Aligned (cSRA) runs must be refused. Uses the local fixture pulled by
/// `pixi run fetch-testdata` (ERR1540848: a small S. pneumoniae cSRA run with a
/// PRIMARY_ALIGNMENT table). Skips if the fixture has not been fetched.
#[test]
fn aligned_run_croaks() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!("skipping aligned_run_croaks: {f} not present (run: pixi run fetch-testdata)");
        return;
    }
    Assert::main_binary()
        .with_args(&[f.as_str()])
        .fails()
        .stderr()
        .contains("aligned")
        .unwrap();
}

/// With --allow-aligned, the same cSRA run is extracted (READ reconstructed from
/// the alignment table) rather than refused. Uses the fetched fixture; skips if
/// absent.
#[test]
fn allow_aligned_extracts() {
    let f = format!(
        "{}/tests/data/ERR1540848/ERR1540848.sra",
        env!("CARGO_MANIFEST_DIR")
    );
    if !std::path::Path::new(&f).exists() {
        eprintln!(
            "skipping allow_aligned_extracts: {f} not present (run: pixi run fetch-testdata)"
        );
        return;
    }
    Assert::main_binary()
        .with_args(&["--allow-aligned", "--single-out", "/dev/null", f.as_str()])
        .succeeds()
        .stdout()
        .contains(">ERR1540848.")
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
