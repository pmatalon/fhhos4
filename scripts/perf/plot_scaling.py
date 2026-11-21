#!/usr/bin/env python3
"""Plots of the strong and weak scaling of U-AMG, setup and solve apart, from the summaries of scaling_mpi.sh: every
summary.txt below the results directories given (those of scripts/cluster/run.pbs, or <dir>/<case>/summary.txt).

Usage: plot_scaling.py <output dir> [<label>=]<results dir>...
Several directories overlay their runs (two commits, the laptop and the cluster...): the last one in solid lines, the
earlier ones dashed, then dotted. Writes strong_<case>.png (top: the cores as OpenMP threads of one process or as MPI
processes of one thread; bottom: all the cores split into MPI processes x OpenMP threads) and weak.png. Needs
matplotlib.
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

LINE = re.compile(r"(\d+) processes x (\d+) threads: setup ([\d.]+) s, solve ([\d.]+) s, (?:(\d+) iterations)?")
THREADS, PROCESSES, SPLIT, IDEAL = "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
STYLES = ("-", "--", ":")  # the last run first
PHASES = (("Setup", 2), ("Solve", 3))

plt.rcParams.update({
    "font.size": 9, "text.color": INK, "axes.labelcolor": INK_2, "axes.edgecolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
})


def parse_summary(path):
    """(processes, threads, setup, solve, iterations) of the runs that ended (a failed run: 999999 s)."""
    rows = []
    for line in path.read_text().splitlines():
        m = LINE.match(line)
        if m and float(m[3]) < 999999:
            rows.append((int(m[1]), int(m[2]), float(m[3]), float(m[4]), int(m[5]) if m[5] else None))
    return rows


def describe(case_dir):
    """'k=0 n=64, 3.64M unknowns', the unknowns from a log of one process (on several, rank 0 prints its own)."""
    m = re.fullmatch(r"k(\d+)n(\d+)", case_dir.name)
    title = f"k={m[1]} n={m[2]}" if m else case_dir.name
    for log in sorted(case_dir.glob("np1_*.log")):
        text = log.read_text(errors="replace")
        mesh = re.search(r"/([a-z0-9]+)_([a-z]+)_N\d+\.msh", text)
        unknowns = re.search(r"Total unknowns: (\d+)", text)
        if mesh:
            title = f"{mesh[1]} {mesh[2]}, {title}"
        if unknowns:
            return f"{title}, {int(unknowns[1]) / 1e6:.2f}M unknowns"
    return title


def load(directory):
    strong, weak = {}, []
    for summary in sorted(directory.rglob("summary.txt")):
        case_dir = summary.parent
        step = case_dir.parent.name
        rows = parse_summary(summary)
        if step == "weak":
            n = int(re.search(r"_n(\d+)", case_dir.name)[1])
            weak += [(p, n, setup, solve, its) for p, t, setup, solve, its in rows]
        elif step in ("native", "sockets"):
            continue  # shortcut: another build or placement of the same cases, not plotted; give them their own runs when needed
        else:
            case = strong.setdefault(case_dir.name, {"rows": [], "title": describe(case_dir)})
            case["rows"] += rows
    return strong, sorted(weak)


def log2_axis(ax, ticks, labels=None):
    ax.set_xscale("log", base=2)
    ax.set_xticks(ticks, labels or [str(t) for t in ticks])
    ax.minorticks_off()


def log_y(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())


def end_labels(ax, labels):
    """Labels (points, text, color) at the end of their lines, moved apart when their ends are close."""
    labels = sorted((points[-1], text, color) for points, text, color in labels if points)
    labels.sort(key=lambda label: label[0][1])
    for j, ((x, y), text, color) in enumerate(labels):
        near = [other[0][1] for other in labels if other[0][0] == x and other[0][1] != y and max(y, other[0][1]) /
                min(y, other[0][1]) < 1.35]
        shift = 0 if not near else (-7 if y < near[0] else 7)
        ax.annotate(text, (x, y), xytext=(6, shift), textcoords="offset points", va="center", fontsize=8,
                    color=INK_2, bbox={"boxstyle": "round,pad=0.15", "fc": "#fcfcfb", "ec": color, "lw": 0.8})


def plot_line(ax, points, color, style, label, latest):
    if points:
        x, y = zip(*points)
        ax.plot(x, y, style, color=color, lw=1.6, marker="o", ms=5, mew=1, label=label,
                mfc=color if latest else "#fcfcfb", mec="#fcfcfb" if latest else color)


def plot_strong(runs, case, out):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.6), constrained_layout=True)
    title = next(r[1][case]["title"] for r in reversed(runs) if case in r[1])
    fig.suptitle(f"Strong scaling: {title}", fontsize=11, x=0.01, ha="left")
    with_case = [r for r in runs if case in r[1]]
    for col, (phase, i) in enumerate(PHASES):
        top, bottom = axes[0][col], axes[1][col]
        cores, splits, totals = set(), [], set()
        for k, (label, strong, _) in enumerate(with_case):
            rows = strong[case]["rows"]
            latest = k == len(with_case) - 1
            style = STYLES[min(len(with_case) - 1 - k, len(STYLES) - 1)]
            suffix = f" ({label})" if len(with_case) > 1 else ""
            threads = sorted((r[1], r[i]) for r in rows if r[0] == 1)
            processes = sorted((r[0], r[i]) for r in rows if r[1] == 1)
            total = max(r[0] * r[1] for r in rows)
            split = sorted((r[0], r[i]) for r in rows if r[0] * r[1] == total)
            cores |= {r[0] * r[1] for r in rows}
            totals.add(total)
            plot_line(top, threads, THREADS, style, "OpenMP" + suffix, latest)
            plot_line(top, processes, PROCESSES, style, "MPI" + suffix, latest)
            plot_line(bottom, split, SPLIT, style, f"{total} cores{suffix}", latest)
            if latest:
                end_labels(top, [(threads, "OpenMP", THREADS), (processes, "MPI", PROCESSES)])
                end_labels(bottom, [(split, f"{total} cores", SPLIT)])
                splits = [(p, f"{p}x{total // p}") for p, _ in split]
                # Ideal: the time of the fewest cores divided by the cores.
                base = min(threads + processes)
                xs = sorted(c for c in cores if c >= base[0])
                top.plot(xs, [base[1] * base[0] / c for c in xs], "--", color=IDEAL, lw=1, label="ideal", zorder=0)
        log2_axis(top, sorted(cores))
        if len(totals) == 1:
            log2_axis(bottom, [p for p, _ in splits], [s for _, s in splits])
            bottom.set_xlabel("MPI processes x OpenMP threads per process")
        else:  # runs on different numbers of cores: the processes alone line up
            log2_axis(bottom, sorted({p for _, strong, _ in with_case for p, t, *_ in strong[case]["rows"]
                                      if p * t == max(r[0] * r[1] for r in strong[case]["rows"])}))
            bottom.set_xlabel("MPI processes (each with cores / processes OpenMP threads)")
        for ax in (top, bottom):
            log_y(ax)
            ax.set_ylabel(f"{phase} (s)")
        top.set_title(f"{phase}: OpenMP (1 process) or MPI (1 thread each)", loc="left", fontsize=10)
        top.set_xlabel("cores")
        bottom.set_title(f"{phase}: all the cores, MPI x OpenMP", loc="left", fontsize=10)
        top.legend(fontsize=8)
        bottom.legend(fontsize=8)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_weak(runs, out):
    with_weak = [r for r in runs if r[2]]
    if not with_weak:
        return False
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.9), constrained_layout=True)
    fig.suptitle("Weak scaling: the same cells per MPI process (n ~ P^(1/3)), one OpenMP thread each", fontsize=11,
                 x=0.01, ha="left")
    for col, (phase, i) in enumerate(PHASES + (("Iterations", 4),)):
        ax = axes[col]
        for k, (label, _, weak) in enumerate(with_weak):
            latest = k == len(with_weak) - 1
            style = STYLES[min(len(with_weak) - 1 - k, len(STYLES) - 1)]
            points = [(r[0], r[i]) for r in weak if r[i] is not None]
            plot_line(ax, points, PROCESSES, style, label, latest)
            if latest:
                if i != 4:
                    ax.axhline(points[0][1], ls="--", color=IDEAL, lw=1, label="ideal", zorder=0)
                log2_axis(ax, [r[0] for r in weak], [f"{r[0]}\nn={r[1]}" for r in weak])
        ax.set_ylim(0, 1.12 * max(r[i] for _, _, weak in with_weak for r in weak if r[i] is not None))
        ax.set_title(phase, loc="left", fontsize=10)
        ax.set_xlabel("MPI processes (mesh size n)")
        ax.set_ylabel(f"{phase} (s)" if i != 4 else phase)
        if col == 0:
            ax.legend(fontsize=8)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return True


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    runs = []
    for arg in sys.argv[2:]:
        label, _, directory = arg.rpartition("=")
        directory = Path(directory)
        if not directory.is_dir():
            sys.exit(f"{directory}: not a directory")
        runs.append((label or directory.name, *load(directory)))
    cases = sorted({case for _, strong, _ in runs for case in strong})
    for case in cases:
        plot_strong(runs, case, out / f"strong_{case}.png")
        print(out / f"strong_{case}.png")
    if plot_weak(runs, out / "weak.png"):
        print(out / "weak.png")


if __name__ == "__main__":
    main()
