#!/usr/bin/env python3
"""
assembly_overview_plot.py

Assembly overview figure for the F. pinicola reference genome.
Three panels: contig structure with rDNA, coverage, heterozygosity.

All input data is read from pipeline output files - no hardcoded values.
Heterozygosity counts come from contig_het_density.tsv (produced by
het_density.py), which uses the permissive definition of heterozygous
sites (all PASS variants where the two genotype alleles differ, including
multiallelic 1/2 calls).

Usage:

This assumes the script is placed in the project root directory.

python3 assembly_overview_plot.py \
    --coverage results/qc/mapping/fpindikaryon/fpindikaryon.coverage.tsv \
    --het_tsv results/qc/het_density_hifiasm/fpindikaryon/contig_het_density.tsv \
    --barrnap results/qc/barrnap/fpindikaryon/barrnap.gff \
    --rrna_contigs results/qc/barrnap/fpindikaryon/rRNA_contigs.txt \
    --output assembly_overview.png

"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import defaultdict
import argparse
import os


def parse_coverage(coverage_path):
    """
    Parse samtools coverage output.
    Columns: rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
    Returns dict: {contig: {"length": int, "coverage": float}}
    """
    contigs = {}
    with open(coverage_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 7:
                continue
            name = cols[0]
            length = int(cols[2])  # endpos
            mean_depth = float(cols[6])
            contigs[name] = {"length": length, "coverage": mean_depth}
    return contigs


def parse_het_counts(tsv_path):
    """
    Read contig_het_density.tsv from het_density.py output.
    Columns: display_label, contig, length_bp, het_count, het_per_Mb
    Returns dict: {contig: count}
    """
    counts = {}
    with open(tsv_path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            contig = parts[1]
            het_count = int(parts[3])
            counts[contig] = het_count
    return counts


def parse_rrna_contigs(rrna_contigs_path):
    """Parse list of rRNA-containing contig names."""
    rrna = set()
    with open(rrna_contigs_path) as f:
        for line in f:
            line = line.strip()
            if line:
                rrna.add(line)
    return rrna


def parse_barrnap(gff_path):
    """Parse barrnap GFF to get rDNA positions per contig."""
    rdna_positions = defaultdict(list)
    if not os.path.exists(gff_path):
        print(f"Warning: barrnap file not found: {gff_path}")
        return rdna_positions

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 5 and cols[2] == "rRNA":
                contig = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                rdna_positions[contig].append((start, end))

    merged = {}
    for contig, intervals in rdna_positions.items():
        intervals.sort()
        result = [intervals[0]]
        for start, end in intervals[1:]:
            if start <= result[-1][1] + 1000:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[contig] = result

    return merged


def classify_contig(length, coverage):
    """Classify contigs by size and coverage."""
    if length >= 1_000_000:
        return "large"
    elif coverage <= 5:
        return "small_low_cov"
    else:
        return "small"


def build_contig_data(coverage_data, het_counts, rrna_set):
    """Combine all sources into per-contig data structure."""
    data = {}
    for name, cov_info in coverage_data.items():
        data[name] = {
            "length": cov_info["length"],
            "coverage": cov_info["coverage"],
            "het_count": het_counts.get(name, 0),
            "rRNA": name in rrna_set,
        }
    return data


def create_overview_plot(contig_data, rdna_positions, output_path):
    sorted_contigs = sorted(
        contig_data.items(), key=lambda x: x[1]["length"], reverse=True
    )
    n_contigs = len(sorted_contigs)

    categories = {}
    for name, data in sorted_contigs:
        categories[name] = classify_contig(data["length"], data["coverage"])

    # Colors
    COL_LARGE = "#4A90B8"
    COL_SMALL = "#D4A373"
    COL_LOW_COV = "#AAAAAA"
    COL_RDNA = "#C0392B"

    bar_colors = {
        "large": COL_LARGE,
        "small": COL_SMALL,
        "small_low_cov": COL_LOW_COV,
    }

    fig, axes = plt.subplots(
        1, 3, figsize=(18, 11),
        gridspec_kw={"width_ratios": [5.5, 1.5, 1.5], "wspace": 0.15}
    )
    ax_main, ax_cov, ax_het = axes

    y_positions = list(range(n_contigs - 1, -1, -1))
    bar_height = 0.65
    max_length = max(d["length"] for _, d in sorted_contigs)

    sep_y = None
    for i, (name, data) in enumerate(sorted_contigs):
        if data["length"] < 1_000_000 and sep_y is None:
            sep_y = y_positions[i] + 0.5
            break

    # Panel A: Contig bars with rDNA
    for i, (name, data) in enumerate(sorted_contigs):
        y = y_positions[i]
        cat = categories[name]
        length = data["length"]
        color = bar_colors[cat]

        ax_main.barh(
            y, length / 1e6, height=bar_height,
            color=color, edgecolor="#333333", linewidth=0.4, alpha=0.9
        )

        if name in rdna_positions:
            for start, end in rdna_positions[name]:
                rdna_x = start / 1e6
                rdna_w = max((end - start) / 1e6, 0.03)
                ax_main.barh(
                    y, rdna_w, left=rdna_x, height=bar_height,
                    color=COL_RDNA, edgecolor="none", alpha=0.9
                )

        suffix = " *" if name.endswith("c") else ""
        ax_main.text(
            -0.12, y, f"{name}{suffix}", ha="right", va="center",
            fontsize=8, fontfamily="monospace"
        )

    if sep_y is not None:
        for ax in axes:
            ax.axhline(sep_y, color="#666666", linestyle="--",
                      linewidth=0.8, alpha=0.6)

    ax_main.set_xlim(0, max_length / 1e6 * 1.05)
    ax_main.set_ylim(-0.8, n_contigs - 0.2)
    ax_main.set_xlabel("Contig length (Mb)", fontsize=12)
    ax_main.set_yticks([])
    ax_main.set_title("Contig structure with rDNA positions",
                      fontsize=13, fontweight="bold", pad=12)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)
    ax_main.spines["left"].set_visible(False)

    # Panel B: Coverage
    cov_values = [d["coverage"] for _, d in sorted_contigs]
    cov_colors_list = [bar_colors[categories[name]] for name, _ in sorted_contigs]

    ax_cov.barh(
        y_positions, cov_values, height=bar_height,
        color=cov_colors_list, edgecolor="#333333", linewidth=0.4, alpha=0.9
    )

    for y, val in zip(y_positions, cov_values):
        ax_cov.text(val + 2, y, f"{val:.0f}×", ha="left", va="center", fontsize=7.5)

    core_covs = [d["coverage"] for n, d in sorted_contigs if categories[n] == "large"]
    median_cov = np.median(core_covs) if core_covs else 0
    if core_covs:
        ax_cov.axvline(median_cov, color="#333333", linestyle=":",
                      linewidth=1.0, alpha=0.5)

    ax_cov.set_xlim(0, max(cov_values) * 1.35)
    ax_cov.set_ylim(-0.8, n_contigs - 0.2)
    ax_cov.set_xlabel("Coverage", fontsize=12)
    ax_cov.set_yticks([])
    ax_cov.set_title("Read depth", fontsize=13, fontweight="bold", pad=12)
    ax_cov.spines["top"].set_visible(False)
    ax_cov.spines["right"].set_visible(False)
    ax_cov.spines["left"].set_visible(False)

    # Panel C: Heterozygosity
    het_values = [
        d["het_count"] / (d["length"] / 1e6)
        for _, d in sorted_contigs
    ]
    het_colors_list = [bar_colors[categories[name]] for name, _ in sorted_contigs]

    ax_het.barh(
        y_positions, het_values, height=bar_height,
        color=het_colors_list, edgecolor="#333333", linewidth=0.4, alpha=0.9
    )

    ax_het.set_xlim(0, max(het_values) * 1.1)
    ax_het.set_ylim(-0.8, n_contigs - 0.2)
    ax_het.set_xlabel("Het. variants per Mb", fontsize=12)
    ax_het.set_yticks([])
    ax_het.set_title("Heterozygosity density", fontsize=13,
                     fontweight="bold", pad=12)
    ax_het.spines["top"].set_visible(False)
    ax_het.spines["right"].set_visible(False)
    ax_het.spines["left"].set_visible(False)

    # Legend
    legend_elements = [
        mpatches.Patch(color=COL_LARGE, label="Contigs ≥ 1 Mb"),
        mpatches.Patch(color=COL_SMALL, label="Contigs < 1 Mb"),
        mpatches.Patch(color=COL_LOW_COV, label="Contigs < 1 Mb, low coverage (≤ 5×)"),
        mpatches.Patch(color=COL_RDNA, label="rDNA regions (barrnap)"),
    ]
    fig.legend(
        handles=legend_elements, loc="lower center", ncol=4,
        fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.01)
    )

    # Summary statistics
    n_large = sum(1 for c in categories.values() if c == "large")
    n_small = sum(1 for c in categories.values() if c in ("small", "small_low_cov"))
    large_size = sum(
        d["length"] for n, d in sorted_contigs if categories[n] == "large"
    ) / 1e6
    total_size = sum(d["length"] for _, d in sorted_contigs) / 1e6
    n_low_cov = sum(1 for c in categories.values() if c == "small_low_cov")
    total_het = sum(d["het_count"] for _, d in sorted_contigs)
    genome_het_per_mb = total_het / total_size

    summary = (
        f"Total: {n_contigs} contigs, {total_size:.1f} Mb   |   "
        f"{n_large} contigs ≥ 1 Mb ({large_size:.1f} Mb)   |   "
        f"{n_small} contigs < 1 Mb (of which {n_low_cov} with ≤ 5× coverage)   |   "
        f"{total_het} het variants ({genome_het_per_mb:.1f}/Mb)"
    )
    annotations_line = (
        f"Dashed line: 1 Mb size threshold   |   "
        f"Dotted line: median read depth of contigs ≥ 1 Mb ({median_cov:.0f}×)"
    )
    fig.suptitle(
        "Assembly structure of the F. pinicola reference genome",
        fontsize=15, fontweight="bold", y=0.99
    )
    fig.text(0.5, 0.955, summary, ha="center", fontsize=10, style="italic")
    fig.text(0.5, 0.93, annotations_line, ha="center", fontsize=9, color="#555555")

    fig.text(
        0.08, -0.015,
        "* circular contig (hifiasm designation)",
        fontsize=8, color="#666666", style="italic"
    )

    plt.tight_layout(rect=[0.08, 0.03, 1, 0.92])
    plt.savefig(output_path, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {output_path}")

    # Print summary table
    print(f"\n{'='*95}")
    print(f"{'Contig':<15} {'Length':>10} {'Category':<14} {'Cov':>6} "
          f"{'Het':>5} {'Het/Mb':>8} {'rDNA%':>7} {'rDNA bp':>10}")
    print(f"{'-'*95}")
    for name, data in sorted_contigs:
        cat = categories[name]
        cov = data["coverage"]
        het_count = data["het_count"]
        het = het_count / (data["length"] / 1e6)
        rdna_bp = 0
        if name in rdna_positions:
            for s, e in rdna_positions[name]:
                rdna_bp += e - s
        rdna_pct = rdna_bp / data["length"] * 100

        print(f"{name:<15} {data['length']:>10,} {cat:<14} {cov:>5.0f}× "
              f"{het_count:>5} {het:>8.1f} {rdna_pct:>6.1f}% {rdna_bp:>10,}")
    print(f"{'-'*95}")
    print(f"{'TOTAL':<15} {int(total_size*1e6):>10,} {'':<14} {'':>6} "
          f"{total_het:>5} {genome_het_per_mb:>8.1f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", required=True,
                        help="samtools coverage output (TSV)")
    parser.add_argument("--het_tsv", required=True,
                        help="contig_het_density.tsv from het_density.py")
    parser.add_argument("--barrnap", required=True,
                        help="barrnap GFF file")
    parser.add_argument("--rrna_contigs", required=True,
                        help="List of rRNA-containing contig names")
    parser.add_argument("--output", default="assembly_overview.png")
    args = parser.parse_args()

    coverage_data = parse_coverage(args.coverage)
    het_counts = parse_het_counts(args.het_tsv)
    rrna_set = parse_rrna_contigs(args.rrna_contigs)
    rdna_positions = parse_barrnap(args.barrnap)

    contig_data = build_contig_data(coverage_data, het_counts, rrna_set)

    create_overview_plot(contig_data, rdna_positions, args.output)
