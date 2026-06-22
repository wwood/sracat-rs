#!/usr/bin/env bash
# Benchmark sracat-rs against `fasterq-dump` (primary comparison) and the
# original C++ `sracat` (reference) on two unaligned runs and an aligned (cSRA)
# run. Inputs are fetched by Snakefile (`pixi run -e download download`).
# Intended to run on a compute node via mqsub, inside the benchmarking pixi env
# (so $CONDA_PREFIX has the ncbi-vdb 2.11 + NGS libs the C++ sracat needs, plus
# fasterq-dump). Each input is first copied to node-local temp ($TMP) so that
# shared-filesystem IO is not the bottleneck (and removed before the next one),
# so size $TMP / memory accordingly:
#
#   mqsub --no-email -t 16 -m 96 --hours 2 -- \
#     pixi run --manifest-path benchmarking/pixi.toml bench
#
# Each tool is timed REPS times (default 4); the first run is a warm-up that the
# plot discards. Writes a machine-readable results.tsv (file, tool, threads, rep,
# seconds, reads, rc) plus a human-readable log to bench_big.txt. Render with
# `pixi run plot` (bars = mean, error bars = stdev over the kept reps).
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"

# Binaries.
RSBIN="${RSBIN:-$REPO/target/release/sracat-rs}"   # built by `pixi run build`
CPPBIN="${CPPBIN:-$HERE/sracat-cpp}"               # built by `pixi run build-cpp`

# Inputs (downloaded by Snakefile into data/<acc>/<acc>.sra). Override via env.
DATA="${DATA:-$HERE/data}"
BIG="${BIG:-$DATA/SRR24704796/SRR24704796.sra}"     # unaligned, ~2.7 GB
MED="${MED:-$DATA/ERR12726217/ERR12726217.sra}"     # unaligned, ~0.6 GB
CSRA="${CSRA:-$DATA/ERR1540848/ERR1540848.sra}"     # aligned cSRA, ~9 MB
T="${T:-16}"        # multi-thread count (each tool is also run at 1 thread)
TIMEOUT="${TIMEOUT:-600}"   # per-command wall limit (s); guards against hangs
REPS="${REPS:-4}"   # runs per tool; the first (rep 1) is a warm-up, dropped in analysis

RESULTS="$HERE/results.tsv"
LIBP="$CONDA_PREFIX/lib64:$CONDA_PREFIX/lib"
# Node-local scratch: prefer PBS's per-job $TMPDIR, else $BENCH_TMP, else /tmp.
TMP="$(mktemp -d "${BENCH_TMP:-${TMPDIR:-/tmp}}/sracatbench.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

exec > >(tee "$HERE/bench_big.txt") 2>&1

printf 'file\ttool\tthreads\trep\tseconds\treads\trc\n' > "$RESULTS"

# run <file_label> <tool> <threads> <count_pat> -- <command...>
# Runs the command REPS times (one TSV row per rep: rep 1 is the warm-up). Stops
# early if a rep fails / times out, so a hung tool isn't repeated REPS times.
run() {
    local flabel="$1" tool="$2" thr="$3" pat="$4"; shift 4
    [ "$1" = "--" ] && shift
    local rep s e rc n secs
    for rep in $(seq 1 "$REPS"); do
        s=$(date +%s.%N)
        timeout "$TIMEOUT" "$@" 2>"$TMP/err" | grep -c "$pat" > "$TMP/cnt"
        rc=${PIPESTATUS[0]}
        e=$(date +%s.%N)
        n=$(cat "$TMP/cnt")
        secs=$(awk -v a="$s" -v b="$e" 'BEGIN{printf "%.2f", b-a}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$flabel" "$tool" "$thr" "$rep" "$secs" "$n" "$rc" >> "$RESULTS"
        printf '%-31s %-22s t=%-2s rep=%d %9ss  rc=%d  reads=%s\n' "$flabel" "$tool" "$thr" "$rep" "$secs" "$rc" "$n"
        if [ "$rc" -ne 0 ]; then
            sed 's/^/    /' "$TMP/err" | tail -3
            break   # don't repeat a failing / timed-out command
        fi
    done
}

echo "host=$(hostname)  T=$T  TMP=$TMP"
echo "RSBIN=$RSBIN"
echo "CPPBIN=$CPPBIN"
echo

# stage <file> : copy <file> (dereferencing symlinks) plus any non-.sra siblings
# (e.g. a cSRA reference) into local temp $TMP/stage, and echo the staged path.
# Reads during the benchmark then hit local storage instead of the shared FS.
# The copy happens before timing starts; caller removes $TMP/stage afterwards.
stage() {
    local f="$1" sib
    rm -rf "$TMP/stage"
    mkdir -p "$TMP/stage"
    cp -L "$f" "$TMP/stage/"
    for sib in "$(dirname "$f")"/*; do
        case "$sib" in
            *.sra) ;;
            *) [ -f "$sib" ] && cp -L "$sib" "$TMP/stage/" ;;
        esac
    done
    echo "$TMP/stage/$(basename "$f")"
}

# Each tool is run at 1 thread and at T threads so the results split cleanly
# into a single-threaded and a multi-threaded comparison. fasterq-dump is the
# primary point of comparison; the C++ sracat is included for reference.

# bench_unaligned <label> <file>
# Unaligned runs: fasterq-dump's sorted mode would build multi-GB temp files, so
# only its (streaming) unsorted mode is compared here.
bench_unaligned() {
    local label="$1" f="$2" tag lf
    if [ ! -r "$f" ]; then
        echo "WARN: input not readable, skipping: $f"
        return
    fi
    tag=$(echo "$label" | tr -c 'A-Za-z0-9' _)
    echo "== $label: $f ($(du -hL "$f" | cut -f1)); staging to $TMP/stage =="
    lf=$(stage "$f")
    run "$label" "sracat-cpp"            1    '^>' -- env LD_LIBRARY_PATH="$LIBP" "$CPPBIN" "$lf"
    run "$label" "fasterq-dump-unsorted" 1    '^>' -- fasterq-dump --fasta-unsorted --stdout -e 1   -t "$TMP/$tag.u1" "$lf"
    run "$label" "fasterq-dump-unsorted" "$T" '^>' -- fasterq-dump --fasta-unsorted --stdout -e "$T" -t "$TMP/$tag.uN" "$lf"
    run "$label" "sracat-rs"             1    '^>' -- env -u LD_LIBRARY_PATH "$RSBIN" --single-out "$TMP/$tag.rs1" "$lf"
    run "$label" "sracat-rs"             "$T" '^>' -- env -u LD_LIBRARY_PATH "$RSBIN" -t "$T" --single-out "$TMP/$tag.rsN" "$lf"
    rm -rf "$TMP/stage"
    echo
}

bench_unaligned "SRR24704796 (2.7 GB, unaligned)" "$BIG"
bench_unaligned "ERR12726217 (0.6 GB, unaligned)" "$MED"

# --------------------------------------------------------------------- cSRA --
# Small aligned run. sracat-rs reconstructs READ per spot with no temp dir
# (default behaviour); fasterq-dump sorted is the order-stable competitor (it
# dumps and sorts via its temp dir), unsorted the streaming one.
if [ -r "$CSRA" ]; then
    L="ERR1540848 (9 MB, cSRA)"
    echo "== $L: $CSRA ($(du -hL "$CSRA" | cut -f1)); staging to $TMP/stage =="
    lf=$(stage "$CSRA")   # also copies the FM211187.1 reference sibling
    run "$L" "sracat-cpp"             1    '^>' -- env LD_LIBRARY_PATH="$LIBP" "$CPPBIN" "$lf"
    run "$L" "fasterq-dump-sorted"    1    '^>' -- fasterq-dump --fasta          --stdout -e 1   -t "$TMP/cfs1" "$lf"
    run "$L" "fasterq-dump-sorted"    "$T" '^>' -- fasterq-dump --fasta          --stdout -e "$T" -t "$TMP/cfsN" "$lf"
    run "$L" "fasterq-dump-unsorted"  1    '^>' -- fasterq-dump --fasta-unsorted --stdout -e 1   -t "$TMP/cfu1" "$lf"
    run "$L" "fasterq-dump-unsorted"  "$T" '^>' -- fasterq-dump --fasta-unsorted --stdout -e "$T" -t "$TMP/cfuN" "$lf"
    run "$L" "sracat-rs"              1    '^>' -- env -u LD_LIBRARY_PATH "$RSBIN" --single-out "$TMP/cs1" "$lf"
    run "$L" "sracat-rs"              "$T" '^>' -- env -u LD_LIBRARY_PATH "$RSBIN" -t "$T" --single-out "$TMP/csN" "$lf"
    rm -rf "$TMP/stage"
    echo
else
    echo "WARN: cSRA input not readable, skipping: $CSRA"
fi

echo "results -> $RESULTS"
echo "DONE"
