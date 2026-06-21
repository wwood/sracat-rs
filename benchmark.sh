#!/usr/bin/env bash
# Benchmark sracat-rs against the C++ sracat and fasterq-dump on one SRA file.
#
#   ./benchmark.sh <file.sra> [threads]
#
# All tools write FASTA to /dev/null; wall time and emitted read counts are
# reported. fasterq-dump is run both unsorted (streaming, multi-threaded, NOT
# order-stable) and sorted (deterministic, builds temp files). sracat-rs is run
# single-threaded (deterministic, no temp) and multi-threaded (deterministic,
# temp files). Large runs should be invoked via mqsub.
set -uo pipefail

SRA="${1:?usage: benchmark.sh <file.sra> [threads]}"
THREADS="${2:-$(nproc)}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
CPP_TOML="$REPO/pixi.toml"
RS_TOML="$HERE/pixi.toml"
RSBIN="$HERE/target/release/sracat-rs"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# run <label> <count-file-or-->  <command...>
# Times the command (stdout -> count file or /dev/null), prints wall seconds.
run() {
    local label="$1" outfile="$2"; shift 2
    local start end
    start=$(date +%s.%N)
    if [ "$outfile" = "-" ]; then
        "$@" >/dev/null 2>"$TMP/err"
    else
        "$@" >"$outfile" 2>"$TMP/err"
    fi
    local rc=$?
    end=$(date +%s.%N)
    local secs; secs=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.2f", b-a}')
    local n=""
    [ "$outfile" != "-" ] && n="reads=$(grep -c '^>' "$outfile" 2>/dev/null)"
    printf "%-42s %8ss  rc=%d  %s\n" "$label" "$secs" "$rc" "$n"
    [ "$rc" -ne 0 ] && sed 's/^/    /' "$TMP/err" | head -3
}

cpp()  { pixi run --manifest-path "$CPP_TOML" -- "$@"; }
fqd()  { pixi run --manifest-path "$CPP_TOML" -- "$@"; }   # fasterq-dump lives in the same env

echo "file    : $SRA ($(du -h "$SRA" 2>/dev/null | cut -f1))"
echo "threads : $THREADS"
echo

run "sracat (C++, 1 thread, all fragments)" "$TMP/cpp.fa"  cpp sracat "$SRA"
run "fasterq-dump --fasta-unsorted (-e1)"        "$TMP/fqu1.fa" \
    fqd fasterq-dump --fasta-unsorted --stdout -e 1 -t "$TMP/fqu1" "$SRA"
run "fasterq-dump --fasta sorted (-e1)"          "$TMP/fqs1.fa" \
    fqd fasterq-dump --fasta --stdout -e 1 -t "$TMP/fqs1" "$SRA"
run "fasterq-dump --fasta-unsorted (-e$THREADS)" "$TMP/fqu.fa" \
    fqd fasterq-dump --fasta-unsorted --stdout -e "$THREADS" -t "$TMP/fqu" "$SRA"
run "fasterq-dump --fasta sorted (-e$THREADS)"   "$TMP/fqs.fa" \
    fqd fasterq-dump --fasta --stdout -e "$THREADS" -t "$TMP/fqs" "$SRA"
run "sracat-rs (1 thread)"               "$TMP/rs1.fa" \
    "$RSBIN" --single-out "$TMP/rs1.single" "$SRA"
run "sracat-rs (-t$THREADS)"             "$TMP/rsn.fa" \
    "$RSBIN" -t "$THREADS" --temp "$TMP" --single-out "$TMP/rsn.single" "$SRA"
