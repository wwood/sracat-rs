#!/usr/bin/env bash
# Big-file benchmark, intended to run on a compute node via mqsub. Writes a
# results table to bench_big.txt in this directory (shared filesystem).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec > >(tee "$HERE/bench_big.txt") 2>&1

CPP=/mnt/weka/pkg/cmr/woodcrob/pixi/cache/envs/sracat-3539029493336913772/envs/default
RSBIN="$HERE/target/release/sracat-rs"
SRACATCPP="$(dirname "$HERE")/sracat"
F=/mnt/weka/scratch/microbiome/woodcrob/non_sensitive/tmp/sylph_staging/SRR24704796/SRR24704796.sra
T=16
TMP=$(mktemp -d)

tm() {
    local s e r
    s=$(date +%s.%N)
    "$@" >/dev/null 2>"$TMP/e"
    r=$?
    e=$(date +%s.%N)
    awk -v a="$s" -v b="$e" -v r="$r" 'BEGIN{printf "%8.1fs rc=%d\n", b-a, r}'
    [ "$r" -ne 0 ] && tail -2 "$TMP/e" | sed 's/^/    /'
}

echo "file: $(du -h "$F" | cut -f1); threads=$T; host=$(hostname)"
echo "== C++ sracat (1thr, all frags) =="
tm env LD_LIBRARY_PATH=$CPP/lib "$SRACATCPP" "$F"
echo "== fasterq-dump --fasta-unsorted -e1 =="
tm env LD_LIBRARY_PATH=$CPP/lib PATH=$CPP/bin:$PATH "$CPP/bin/fasterq-dump" --fasta-unsorted --stdout -e 1 -t "$TMP/fqu1" "$F"
echo "== fasterq-dump --fasta-unsorted -e$T =="
tm env LD_LIBRARY_PATH=$CPP/lib PATH=$CPP/bin:$PATH "$CPP/bin/fasterq-dump" --fasta-unsorted --stdout -e "$T" -t "$TMP/fqu" "$F"
echo "== sracat-rs -t1 =="
tm env -u LD_LIBRARY_PATH "$RSBIN" --single-out "$TMP/s1" "$F"
echo "== sracat-rs -t$T =="
tm env -u LD_LIBRARY_PATH "$RSBIN" -t "$T" --temp "$TMP" --single-out "$TMP/sn" "$F"

rm -rf "$TMP"
echo "DONE"
