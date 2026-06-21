# Pull SRA test fixtures used by the integration tests.
#
# ERR1540848 is a small (~8 MB) Streptococcus pneumoniae cSRA (aligned) run,
# aligned to a ~2 Mb reference (FM211187.1). It is too large to commit, so it is
# fetched on demand here and git-ignored. tests/cli.rs::aligned_run_croaks uses
# it to check that sracat-rs refuses aligned runs.
#
# Run via the pixi env (prefetch comes from sra-tools in pixi.toml):
#   pixi run fetch-testdata


rule all:
    input:
        "tests/data/ERR1540848/ERR1540848.sra",


rule prefetch_aligned_csra:
    output:
        "tests/data/ERR1540848/ERR1540848.sra",
    log:
        "logs/prefetch_ERR1540848.log",
    benchmark:
        "benchmarks/prefetch_ERR1540848.tsv"
    shell:
        "prefetch --output-directory tests/data ERR1540848 > {log} 2>&1"
