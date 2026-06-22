# Pull SRA test fixtures used by the integration tests. Both are too large to
# commit, so they are fetched on demand here and git-ignored.
#
#   ERR1540848 - a small (~9 MB) Streptococcus pneumoniae cSRA (aligned) run,
#                aligned to a ~2 Mb reference (FM211187.1). Exercises the
#                aligned-run handling.
#   DRR033172  - a small (~0.4 MB) unaligned WGS paired run (46,282 genuine
#                read pairs, spanning several decode chunks). Exercises
#                multi-threaded -1/-2 pair splitting.
#
# Run via the pixi env (prefetch comes from sra-tools in pixi.toml):
#   pixi run fetch-testdata


rule all:
    input:
        "tests/data/ERR1540848/ERR1540848.sra",
        "tests/data/DRR033172/DRR033172.sra",


rule prefetch_aligned_csra:
    output:
        "tests/data/ERR1540848/ERR1540848.sra",
    log:
        "logs/prefetch_ERR1540848.log",
    benchmark:
        "benchmarks/prefetch_ERR1540848.tsv"
    shell:
        "prefetch --output-directory tests/data ERR1540848 > {log} 2>&1"


rule prefetch_unaligned_paired:
    output:
        "tests/data/DRR033172/DRR033172.sra",
    log:
        "logs/prefetch_DRR033172.log",
    benchmark:
        "benchmarks/prefetch_DRR033172.tsv"
    shell:
        "prefetch --output-directory tests/data DRR033172 > {log} 2>&1"
