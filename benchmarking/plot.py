#!/usr/bin/env python
"""Render results.tsv (from bigbench.sh) into per-thread-count SVG charts.

Each tool is timed several times per run; the first rep is a warm-up and is
dropped. Bars show the mean wall-clock time over the remaining reps, with error
bars giving the sample standard deviation. One SVG per thread count (e.g.
comparison-t1.svg, comparison-t16.svg); within each, one panel per input file.
fasterq-dump is the primary point of comparison (coloured distinctly), sracat-rs
is highlighted, and the C++ sracat is greyed.
"""

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.tsv"

RS = "#1b7837"  # sracat-rs (highlight)
FQD = "#2c7fb8"  # fasterq-dump (primary comparison)
CPP = "#999999"  # C++ sracat (reference)


def colour(tool: str) -> str:
    if tool == "sracat-rs":
        return RS
    if tool.startswith("fasterq-dump"):
        return FQD
    return CPP


def render(agg: pl.DataFrame, threads: int, files: list[str]) -> Path:
    sub_all = agg.filter(pl.col("threads") == threads)
    panels = [f for f in files if f in sub_all["file"].unique().to_list()]
    fig, axes = plt.subplots(
        len(panels), 1,
        figsize=(9, 1.0 + 1.4 * len(panels) + 0.4 * sub_all.height),
        squeeze=False,
    )
    for ax, fname in zip(axes[:, 0], panels):
        sub = sub_all.filter(pl.col("file") == fname).sort("mean", descending=True)
        names = sub["tool"].to_list()
        means = sub["mean"].to_list()
        stds = sub["std"].to_list()
        colours = [colour(t) for t in names]
        bars = ax.barh(
            names, means, xerr=stds, color=colours,
            capsize=4, error_kw={"ecolor": "#333333", "lw": 1},
        )
        labels = [f"{m:.1f}±{s:.1f}s" for m, s in zip(means, stds)]
        ax.bar_label(bars, labels=labels, padding=3, fontsize=8)
        ax.set_xlabel("wall-clock seconds (lower is better)")
        ax.set_title(f"{fname} run", loc="left", fontweight="bold")
        ax.margins(x=0.20)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=RS),
        plt.Rectangle((0, 0), 1, 1, color=FQD),
        plt.Rectangle((0, 0), 1, 1, color=CPP),
    ]
    fig.suptitle(
        f"sracat-rs vs fasterq-dump vs sracat (C++) — {threads} thread"
        f"{'s' if threads != 1 else ''}",
        fontweight="bold",
    )
    fig.legend(
        handles,
        ["sracat-rs", "fasterq-dump", "sracat (C++)"],
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = HERE / f"comparison-t{threads}.svg"
    fig.savefig(out, format="svg")
    plt.close(fig)
    return out


def main() -> None:
    raw = pl.read_csv(RESULTS, separator="\t")
    # Drop the warm-up (rep 1) and any failed reps, then average per tool/run.
    kept = raw.filter((pl.col("rep") >= 2) & (pl.col("rc") == 0))
    if kept.height == 0:
        raise SystemExit(f"no successful non-warm-up rows in {RESULTS}")
    agg = kept.group_by(["file", "tool", "threads"]).agg(
        pl.col("seconds").mean().alias("mean"),
        pl.col("seconds").std(ddof=1).fill_null(0.0).alias("std"),
        pl.len().alias("n"),
    )
    files = raw["file"].unique(maintain_order=True).to_list()
    for threads in sorted(agg["threads"].unique().to_list()):
        out = render(agg, threads, files)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
