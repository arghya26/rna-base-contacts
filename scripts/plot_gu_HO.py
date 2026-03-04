#!/usr/bin/env python3
"""
Analyze GU NH···O contacts from gu_HO_summary.txt

Produces:
1. Bar chart: frequency vs H–O contact type (consistent colors with lines).
2. Line-style (step) plots: distance distributions per H–O type.

Usage:
    python plot_gu_HO.py gu_HO_summary.txt
"""

import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def read_HO_summary(path):
    """
    Returns:
      counts: Counter of H_atom-O_atom strings
      distances_by_type: dict[type] -> list of distances (float)
    """
    counts = Counter()
    distances_by_type = defaultdict(list)

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 12:
                continue
            H_atom = fields[5].strip()
            O_atom = fields[10].strip()
            dist = float(fields[11])

            ho_type = f"{H_atom}-{O_atom}"
            counts[ho_type] += 1
            distances_by_type[ho_type].append(dist)

    return counts, distances_by_type


def assign_colors(counts):
    """
    Assign a consistent color per H–O type, using a matplotlib colormap.
    Returns: dict[type] -> color
    """
    types = [t for t, _ in counts.most_common()]
    cmap = plt.get_cmap("tab10")  # good qualitative palette
    colors = {}
    for i, t in enumerate(types):
        colors[t] = cmap(i % 10)
    return colors


def plot_type_frequency(counts, colors, out_png="HO_type_frequency.png"):
    # Sort types by descending count
    types = [t for t, _ in counts.most_common()]
    freqs = [counts[t] for t in types]
    bar_colors = [colors[t] for t in types]

    plt.figure(figsize=(8, 4))
    plt.bar(range(len(types)), freqs, color=bar_colors)
    plt.xticks(range(len(types)), types, rotation=45, ha="right")
    plt.ylabel("Count")
    plt.title("Frequency of NH–O contact types in GU pairs")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_distance_bars(distances_by_type, colors, min_count=5,
                        out_png="HO_distance_histograms_stepstyle.png"):
    """
    Transparent bar/step histograms of distance distributions for each H–O type.
    Only types with at least min_count occurrences are plotted.
    """
    plt.figure(figsize=(8, 5))
    ax = plt.gca()

    # Minor ticks every 0.1 Å on x axis
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1)) # MultipleLocator(0.1) enforces minor ticks every 0.1 on the x axis
    ax.tick_params(axis="x", which="minor", length=3)       # optional: style of minor ticks
    
    ### Collect all distances to define a common bin range ###
    all_dists = [d for dists in distances_by_type.values() for d in dists]
    if not all_dists:
        print("No distances found; skipping distance plot.")
        return

    d_min = min(all_dists) - 0.1
    d_max = max(all_dists) + 0.1
    bins = np.linspace(d_min, d_max, 30)
    
    #bin_width = 0.075     # Å
    #bins = np.arange(d_min, d_max + bin_width, bin_width)
    #print("bins:", bins)
    #print("n_bins:", len(bins) - 1)      # number of histogram bins

    for ho_type, dists in distances_by_type.items():
        if len(dists) < min_count:
            continue
        color = colors.get(ho_type, "gray")
        plt.hist(
            dists,
            bins=bins,
            alpha=0.40,                  # transparency like image.jpg
            color=color,
            edgecolor="#999999",
            linewidth=0.30,
            label=f"{ho_type} (n={len(dists)})"
        )

    plt.xlabel("H–O distance (Å)")
    plt.ylabel("Count")
    plt.title("Distance distributions for NH–O contact types")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(out_png, dpi=600)
    plt.close()


def plot_distance_lines(distances_by_type, colors, min_count=5,
                        out_png="HO_distance_histograms_linestyle.png"):
    """
    Line plots of distance distributions for each H–O type.
    Only types with at least min_count occurrences are plotted.
    """
    plt.figure(figsize=(8, 5))
    ax = plt.gca()

    # Minor ticks every 0.1 Å on x axis
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1)) # MultipleLocator(0.1) enforces minor ticks every 0.1 on the x axis
    ax.tick_params(axis="x", which="minor", length=3)       # optional: style of minor ticks

    # Collect all distances to define a common bin range
    all_dists = [d for dists in distances_by_type.values() for d in dists]
    if not all_dists:
        print("No distances found; skipping distance plot.")
        return

    d_min = min(all_dists) - 0.1
    d_max = max(all_dists) + 0.1
    bins = np.linspace(d_min, d_max, 30)

    for ho_type, dists in distances_by_type.items():
        if len(dists) < min_count:
            continue
        color = colors.get(ho_type, "gray")
        hist, edges = np.histogram(dists, bins=bins)
        centers = 0.5 * (edges[:-1] + edges[1:])
        # Line-style curve (like your example)
        plt.plot(centers, hist, "-", color=color, label=f"{ho_type} (n={len(dists)})")

    plt.xlabel("H–O distance (Å)")
    plt.ylabel("Count")
    plt.title("Distance distributions for NH–O contact types")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(out_png, dpi=600)
    plt.close()


def main():
    if len(sys.argv) != 2:
        print("Usage: python analyze_gu_HO.py gu_HO_summary.txt")
        sys.exit(1)

    path = sys.argv[1]
    counts, distances_by_type = read_HO_summary(path)

    if not counts:
        print("No H–O contacts found.")
        sys.exit(0)

    colors = assign_colors(counts)

    # 1. Frequency vs H–O type (bars, colored by type)
    plot_type_frequency(counts, colors)

    # 2. Distance-based bar/hist plots, colored by H–O type
    plot_distance_bars(distances_by_type, colors, min_count=5) 

    # 2. Distance-based line plots, colored by H–O type
    plot_distance_lines(distances_by_type, colors, min_count=5)


if __name__ == "__main__":
    main()
