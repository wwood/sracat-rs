# sracat-rs

A fast, deterministic reimplementation of [`sracat`](../) in Rust, reading SRA
files directly through the **ncbi-vdb cursor C API** (via a small C shim) rather
than the legacy NGS C++ engine.

## Why

The original `sracat` uses the NGS C++ API (`ncbi::NGS::openReadCollection`),
which is only available in ncbi-vdb 2.x. `sracat-rs` instead binds the
lower-level VDB cursor API present in current ncbi-vdb (3.4.1), reading the
`SEQUENCE` table's `READ` / `READ_LEN` / `READ_TYPE` columns in storage order.

Benefits:

- **Repeatable order.** Reads are emitted in `SEQUENCE`-table row order,
  single-threaded — byte-identical across runs (unlike `fasterq-dump
  --fasta-unsorted` with multiple threads).
- **Faster.** ~3–4× quicker than the NGS-based C++ `sracat` on the same data
  (no per-fragment allocations, no NGS virtual-dispatch layer; reads column
  blobs directly).
- **Paired / single aware.** Spots with two biological reads are emitted
  interleaved (`/1`, `/2`); spots with a single biological read are routed to a
  separate stream.
- **Refuses aligned runs.** If the run is aligned (a `PRIMARY_ALIGNMENT` table
  is present), `READ` in `SEQUENCE` is reconstructed from alignments rather than
  stored; `sracat-rs` errors out instead of silently doing expensive/incorrect
  work.

## Build

Everything (Rust toolchain, C compiler, ncbi-vdb headers + shared library) is
provided by the pixi environment:

```sh
pixi run build      # cargo build --release
pixi run test       # cargo test
```

The release binary lands at `target/release/sracat-rs`. The path to
`libncbi-vdb.so` is baked in as an rpath, so the binary runs outside `pixi run`
as long as that environment still exists.

## Usage

```
sracat-rs [OPTIONS] <SRA>...

  -o, --output-prefix <PREFIX>   write pairs to <PREFIX>.paired.{fasta,fastq}
                                 and singles to <PREFIX>.single.{fasta,fastq}
      --single-out <FILE>        when streaming pairs to stdout, write
                                 single/orphan reads here
      --qual                     write FASTQ (with quality) instead of FASTA
      --include-technical        include technical reads (default: biological only)
  -t, --threads <N>              parallel extraction threads (default 1)
      --temp <DIR>               temp dir for --threads > 1
```

Default behaviour streams interleaved paired reads to **stdout**. If the run
contains any unpaired (single/orphan) reads and no destination for them is
given, `sracat-rs` refuses rather than dropping them — pass `--single-out` or
`-o`.

```sh
# stream interleaved pairs to stdout
sracat-rs run.sra | head

# split paired and single output into files
sracat-rs -o out run.sra            # -> out.paired.fasta, out.single.fasta

# a single-end run, streaming pairs (none) to stdout and singles to a file
sracat-rs --single-out singles.fasta run.sra

# FASTQ, split output, 16 threads
sracat-rs --qual -t 16 -o out run.sra
```

## Threads

`--threads N` (N > 1) decodes contiguous row ranges in parallel — each worker
owns a cursor and formats chunks into memory buffers, which a single writer
thread emits in chunk order through a bounded channel (no temp files). The
output is **byte-identical** to the single-threaded run (verified by md5 for
both FASTA and FASTQ). Decoding (especially `QUALITY`) is the CPU bottleneck;
formatting/writing is only ~15–20% of the time, so this parallelises the part
that matters.

Benchmark on a dedicated allocation (not a shared login node — CPU contention
there hides scaling entirely).

## Benchmarks

`./benchmark.sh <file.sra> [threads]` times `sracat-rs` against the C++ `sracat`
and `fasterq-dump`, reporting wall time and read counts.

**2.7 GB run (SRR24704796), 8 dedicated cores:**

| tool                                | time    | speedup |
| ----------------------------------- | ------- | ------- |
| C++ sracat (1 thread)               | 40.9 s  | —       |
| fasterq-dump `--fasta-unsorted` -e1 | 24.2 s  |         |
| fasterq-dump `--fasta-unsorted` -e8 | 7.5 s   |         |
| sracat-rs -t1                       | 21.8 s  | 1.0×    |
| sracat-rs -t2                       | 11.3 s  | 1.9×    |
| sracat-rs -t4                       | 5.8 s   | 3.8×    |
| **sracat-rs -t8**                   | 3.5 s   | 6.2×    |

`sracat-rs` scales near-linearly to the 8 available cores. At -t8 it is ~2.1×
faster than `fasterq-dump -e8` and ~12× faster than the single-threaded C++
`sracat`; even -t1 beats the C++ tool (~1.9×) and matches single-threaded
fasterq-dump — while being deterministic and running against current ncbi-vdb
(3.x). (Benchmark on dedicated cores: on a contended shared node, scaling is
masked by CPU competition.)
