#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import gzip
import math
from collections import defaultdict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def read_fai(fai_path):
    contigs = []
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            length = int(parts[1])
            contigs.append((name, length))
    return contigs


def is_het_gt(gt_str):
    """
    Return True if genotype appears heterozygous.
    Supports both phased (0|1) and unphased (0/1).
    """
    if gt_str is None or gt_str == ".":
        return False
    sep = "|" if "|" in gt_str else ("/" if "/" in gt_str else None)
    if sep is None:
        return False
    a, b = gt_str.split(sep, 1)
    if a == "." or b == ".":
        return False
    return a != b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--out_contigs", required=True)
    ap.add_argument("--out_windows", required=True)
    ap.add_argument("--plot_contigs", required=True)
    ap.add_argument("--plot_windows", required=True)
    args = ap.parse_args()

    contigs = read_fai(args.fai)
    contig_len = {c: L for c, L in contigs}

    het_count = defaultdict(int)
    win_counts = defaultdict(int)

    sample_col = None
    format_col = None

    with open_maybe_gz(args.vcf) as f:
        for line in f:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if len(header) >= 10:
                    format_col = 8
                    sample_col = 9
                else:
                    raise RuntimeError("VCF appears to lack a sample column (no genotype field).")
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            filt = parts[6]
            if filt != "PASS" and filt != ".":
                continue
            
            chrom = parts[0]
            if chrom not in contig_len:
                continue


            pos = int(parts[1])
            fmt = parts[format_col].split(":")
            sample = parts[sample_col].split(":")

            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue

            gt = sample[gt_idx] if gt_idx < len(sample) else None
            if not is_het_gt(gt):
                continue

            het_count[chrom] += 1
            pos0 = pos - 1
            win_idx = pos0 // args.window
            win_counts[(chrom, win_idx)] += 1

    # Build contig summary
    rows = []
    for chrom, L in contigs:
        c = het_count.get(chrom, 0)
        per_mb = c / (L / 1e6) if L > 0 else 0.0
        rows.append((chrom, L, c, per_mb))

    # Sort by heterozygous variants per Mb, then by total count
    rows_sorted = sorted(rows, key=lambda x: (x[3], x[2]), reverse=True)

    # Write contig table with display labels
    with open(args.out_contigs, "w") as out:
        out.write("display_label\tcontig\tlength_bp\thet_count\thet_per_Mb\n")
        for i, (chrom, L, c, per_mb) in enumerate(rows_sorted, start=1):
            out.write(f"C{i}\t{chrom}\t{L}\t{c}\t{per_mb:.3f}\n")

    # Write window table
    with open(args.out_windows, "w") as out:
        out.write("contig\twindow_start\twindow_end\thet_count\n")
        for chrom, L in contigs:
            nwin = math.ceil(L / args.window)
            for w in range(nwin):
                start = w * args.window
                end = min((w + 1) * args.window, L)
                c = win_counts.get((chrom, w), 0)
                out.write(f"{chrom}\t{start}\t{end}\t{c}\n")

    # Plot 1: top contigs by heterozygous variants per Mb
    top_rows = rows_sorted[:20]
    labels = [f"C{i+1}" for i in range(len(top_rows))]
    values = [r[3] for r in top_rows]

    plt.figure(figsize=(12, 6))
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(values)), labels, rotation=0)
    plt.ylabel("Heterozygous variants per Mb")
    plt.xlabel("Contigs ranked by heterozygosity density")
    plt.title("Heterozygosity density across contigs")
    plt.tight_layout()
    plt.savefig(args.plot_contigs, dpi=300)
    plt.close()

    # Plot 2: windows along top-N contigs, still using real contigs but labeled C1..Cn
    top_contigs = [r[0] for r in rows_sorted[:args.top]]
    label_map = {chrom: f"C{i+1}" for i, chrom in enumerate(top_contigs)}

    x = []
    y = []
    tick_pos = []
    tick_lab = []
    offset = 0

    for chrom in top_contigs:
        L = contig_len[chrom]
        nwin = math.ceil(L / args.window)
        midpoints = []
        counts = []

        for w in range(nwin):
            start = w * args.window
            end = min((w + 1) * args.window, L)
            mid = offset + (start + end) / 2
            c = win_counts.get((chrom, w), 0)
            midpoints.append(mid)
            counts.append(c)

        if midpoints:
            x.extend(midpoints)
            y.extend(counts)
            tick_pos.append(offset + L / 2)
            tick_lab.append(label_map[chrom])

        offset += L + args.window

    plt.figure(figsize=(14, 6))
    plt.scatter(x, y, s=8)
    plt.ylabel(f"Heterozygous variants per {args.window // 1000} kb window")
    plt.xlabel("Top contigs ranked by heterozygosity density")
    plt.title("Distribution of heterozygous variants across top contigs")
    plt.xticks(tick_pos, tick_lab, rotation=0)
    plt.tight_layout()
    plt.savefig(args.plot_windows, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
