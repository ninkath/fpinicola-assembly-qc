# F. pinicola assembly, QC and heterozygosity pipeline

Snakemake workflow for de novo assembly, polishing, evaluation and inter-nuclear
heterozygosity analysis of Oxford Nanopore long-read data from a dikaryotic
*Fomitopsis pinicola* isolate.

This workflow was developed for the master's thesis *De novo assembly and
annotation of a near-chromosome-level reference genome for Fomitopsis pinicola,
and comparative genomic analysis with closely related species* (Nina Thorstensen,
University of South-Eastern Norway, 2026).

A companion workflow for genome annotation and comparative analysis is available
[[LINK to annotation/comparative repository here](https://github.com/ninkath/fpinicola-annotation-comparative)]

## What this workflow does

Starting from basecalled Oxford Nanopore FASTQ files, the workflow produces:

- **Filtered read sets** for assembly (Filtlong, ≥ 20 kb, target 2.5 Gb) and for
  variant calling (≥ 5 kb, target 10 Gb)
- **De novo assemblies** with both Flye and hifiasm, for direct comparison
- **Polished consensus** from hifiasm followed by two rounds of Racon and one
  round of Medaka
- **Assembly QC** at each stage (QUAST contiguity, BUSCO completeness, Merqury
  k-mer-based quality value, samtools coverage, mosdepth windowed depth)
- **Telomeric repeat detection** with tidk, including discovery of candidate
  motifs and search for canonical and variant fungal telomere sequences
- **Heterozygosity analysis** with Clair3 on the variant-calling read set,
  including per-contig and 100 kb-windowed density tracks
- **Per-window tracks** used as input for the circos visualisation in the thesis
  (GC content, read coverage, repeat composition, gene density, heterozygosity,
  telomere caps)

## Requirements

- **Snakemake** ≥ 9.13.4
- **Apptainer** (or Singularity) for containerised tools
- **Conda** or Micromamba for the small number of plotting environments
- A local **BUSCO database** for the `polyporales_odb12` lineage, placed under
  `resources/busco_db/` (see `qc.busco_db_path` in `config.yaml`)
- A local **Clair3 model** for ONT R10.4.1 super-accurate basecalls
  (`r1041_e82_400bps_sup_v500` or equivalent). The default `clair3_model` path
  in `config.yaml` points to the model bundled in the Clair3 container; if you
  use a different model, you may need to update this path.

## Repository contents

```
.
├── Snakefile                        # Snakemake workflow definition
├── config.yaml                      # Sample paths, parameters, container images
├── envs/
│   └── hetplot.yaml                 # Conda environment for heterozygosity plotting
├── scripts/
│   ├── assembly_overview_plot.py    # Script for the assembly overview plot.
│   └── het_density.py               # Per-contig and windowed heterozygosity from Clair3 VCF
├── README.md                        # This file
└── LICENSE
```

The `results/` directory is created at runtime and is not version-controlled.
to live outside the repository.

## Configuration

Open `config.yaml` and update the sample path:

```yaml
samples:
  fpindikaryon: "data/reads/fpindikaryon.fastq.gz"
```

Other parameters (Filtlong thresholds, polishing rounds, BUSCO lineage, container
images) are pre-set to the values used in the thesis.

## Running the workflow

This workflow does **not** define a default `rule all`. Each result is generated
by specifying its target output file. Rules are intended to be run in the order
described in the thesis methods.

A typical run sequence:

```bash
# 1. Read QC and filtering
snakemake --use-apptainer --cores 8 \
    results/qc/nanoplot_raw/fpindikaryon/NanoStats.txt \
    results/filtered/fpindikaryon_clean.fastq.gz \
    results/qc/nanoplot_filtered/fpindikaryon/NanoStats.txt

# 2. Assembly with both assemblers
snakemake --use-apptainer --cores 16 \
    results/assembly/flye/fpindikaryon/assembly.fasta \
    results/assembly/hifiasm/fpindikaryon/assembly.fasta

# 3. Unpolished QC
snakemake --use-apptainer --cores 8 \
    results/qc/quast/unpolished_hifiasm/fpindikaryon/report.txt \
    results/qc/busco/unpolished_hifiasm/fpindikaryon/short_summary.txt

# 4. Polishing (Racon + Medaka)
snakemake --use-apptainer --cores 8 \
    results/medaka/fpindikaryon/consensus.fasta

# 5. Polished assembly QC and read-back mapping
snakemake --use-apptainer --cores 8 \
    results/qc/quast_polished/fpindikaryon/report.txt \
    results/qc/busco_polished/fpindikaryon/short_summary.txt \
    results/qc/quast_core/fpindikaryon/report.txt \
    results/qc/mapping/fpindikaryon/fpindikaryon.flagstat.txt \
    results/merqury/fpindikaryon/merqury_out.qv

# 6. Telomere search and per-window tracks
snakemake --use-apptainer --cores 4 \
    results/circos/fpindikaryon/tidk_explore_top.tsv \
    results/circos/fpindikaryon/tidk_ttagg_telomeric_repeat_windows.tsv \
    results/circos/fpindikaryon/coverage.regions.bed.gz \
    results/circos/fpindikaryon/gc_content_10kb.tsv

# 7. Heterozygosity analysis on the variant-calling read set
snakemake --use-apptainer --cores 8 \
    results/filtered_phasing/fpindikaryon_clean.fastq.gz \
    results/clair3_hifiasm/fpindikaryon/merge_output.vcf.gz \
    results/clair3_hifiasm/fpindikaryon/het_stats.txt

snakemake --use-apptainer --use-conda --cores 4 \
    results/qc/het_density_hifiasm/fpindikaryon/contig_het_density.tsv
```

A dry-run of any command above (add `-n`) shows the rules that will execute
without producing files.

## Assembly overview plot script

This script should be run from the project root, and paths for input files should be confirmed before running.

## License

MIT License (see LICENSE).
