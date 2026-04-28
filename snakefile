# ============================================================================
# Snakefile: Assembly, polishing, QC, and heterozygosity pipeline
# ============================================================================
configfile: "config.yaml"

def ref_fasta(w):
    ref_choice = config.get("reference_fasta", {}).get("use", "polished")
    if ref_choice == "polished":
        return config["paths"]["polished_ref"].format(sample=w.sample)
    elif ref_choice == "unpolished":
        return config["paths"]["unpolished_ref"].format(sample=w.sample)
    else:
        raise ValueError(f"Unknown reference_fasta.use: {ref_choice}")

SAMPLES = list(config["samples"].keys())
ASSEMBLERS = ["flye", "hifiasm"]

# No default rule_all target is defined in this workflow.
# Run individual rules manually in the desired order.

# ============================================================================
# READ PREPARATION
# ============================================================================
# nanoplot_raw: read length and quality summary for the raw fastq
rule nanoplot_raw:
    input:
        reads=lambda w: config["samples"][w.sample]
    output:
        report="results/qc/nanoplot_raw/{sample}/NanoPlot-report.html",
        stats="results/qc/nanoplot_raw/{sample}/NanoStats.txt"
    threads: config["threads"]["medium"]
    container: config["containers"]["nanoplot"]
    shell:
        r"""
        set -euo pipefail
        export LC_ALL=C
        mkdir -p $(dirname {output.report})
        NanoPlot --fastq {input.reads} -o $(dirname {output.report}) -t {threads} \
            --huge --tsv_stats --loglength --N50 --plots hex dot
        test -s {output.stats}
        """

# filter_reference: Filtlong filter for assembly (length 20 kb, target 2.5 Gb)
rule filter_reference:
    input:
        fq=lambda w: config["samples"][w.sample]
    output:
        fq="results/filtered/{sample}_clean.fastq.gz"
    params:
        min_length=config["filtering"]["reference"]["min_length"],
        target_bases=config["filtering"]["reference"]["target_bases"]
    container: config["containers"]["filtlong"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fq})

        export LC_ALL=C
        export LANG=C

        filtlong --min_length {params.min_length} \
                 --target_bases {params.target_bases} \
                 {input.fq} | gzip > {output.fq}

        test -s {output.fq}
        """
        
# nanoplot_filtered: stats after assembly filtering
rule nanoplot_filtered:
    input:
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        report="results/qc/nanoplot_filtered/{sample}/NanoPlot-report.html",
        stats="results/qc/nanoplot_filtered/{sample}/NanoStats.txt"
    threads: config["threads"]["medium"]
    container: config["containers"]["nanoplot"]
    shell:
        r"""
        set -euo pipefail
        export LC_ALL=C
        mkdir -p $(dirname {output.report})
        NanoPlot --fastq {input.reads} -o $(dirname {output.report}) -t {threads} \
            --tsv_stats --loglength --N50 --plots hex dot
        test -s {output.stats}
        """

# filter_phasing: less stringent filter for variant calling (length 5 kb, target 10 Gb)
rule filter_phasing:
    input:
        fq=lambda w: config["samples"][w.sample]
    output:
        fq="results/filtered_phasing/{sample}_clean.fastq.gz"
    params:
        min_length=config["filtering"]["phasing"]["min_length"],
        target_bases=config["filtering"]["phasing"]["target_bases"]
    container: config["containers"]["filtlong"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fq})

        export LC_ALL=C
        export LANG=C

        filtlong --min_length {params.min_length} \
                 --target_bases {params.target_bases} \
                 {input.fq} | gzip > {output.fq}

        test -s {output.fq}
        """

# nanoplot_filtered_phasing: stats after variant-calling filter
rule nanoplot_filtered_phasing:
    input:
        reads="results/filtered_phasing/{sample}_clean.fastq.gz"
    output:
        report="results/qc/nanoplot_filtered_phasing/{sample}/NanoPlot-report.html",
        stats="results/qc/nanoplot_filtered_phasing/{sample}/NanoStats.txt"
    threads: config["threads"]["medium"]
    container: config["containers"]["nanoplot"]
    shell:
        r"""
        set -euo pipefail
        export LC_ALL=C
        mkdir -p $(dirname {output.report})
        NanoPlot --fastq {input.reads} -o $(dirname {output.report}) -t {threads} \
            --tsv_stats --loglength --N50 --plots hex dot
        test -s {output.stats}
        """

# ============================================================================
# DE NOVO ASSEMBLY (unpolished)
# ============================================================================
# assemble_flye: de novo assembly with Flye
rule assemble_flye:
    input:
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        fasta="results/assembly/flye/{sample}/assembly.fasta",
        gfa="results/assembly/flye/{sample}/assembly_graph.gfa"
    params:
        genome_size=config["assembly"]["genome_size"],
        mode=config["assembly"]["flye_mode"],
        outdir=lambda w: f"results/assembly/flye/{w.sample}"
    threads: config["assembly"]["threads"]
    container: config["containers"]["flye"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        flye {params.mode} "{input.reads}" \
             --genome-size "{params.genome_size}" \
             --out-dir "{params.outdir}" \
             --threads {threads}
        test -s "{output.fasta}"
        """

# assemble_hifiasm: de novo assembly with hifiasm in ONT mode
rule assemble_hifiasm:
    input:
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        fasta="results/assembly/hifiasm/{sample}/assembly.fasta"
    params:
        outdir=lambda w: f"results/assembly/hifiasm/{w.sample}",
        prefix=lambda w: f"results/assembly/hifiasm/{w.sample}/hifiasm_out"
    threads: config["assembly"]["threads"]
    container: config["containers"]["hifiasm"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        hifiasm -o "{params.prefix}" -t {threads} --ont "{input.reads}"
        awk '/^S/{{print ">"$2"\n"$3}}' "{params.prefix}.bp.p_ctg.gfa" > "{output.fasta}"
        test -s "{output.fasta}"
        """

# ============================================================================
# UNPOLISHED QC (BUSCO + QUAST)
# ============================================================================
# busco_unpolished: BUSCO completeness on raw unpolished assembly
rule busco_unpolished:
    input:
        assembly="results/assembly/{assembler}/{sample}/assembly.fasta"
    output:
        summary="results/qc/busco/unpolished_{assembler}/{sample}/short_summary.txt"
    params:
        lineage=config["qc"]["busco_lineage"],
        db_path=config["qc"]["busco_db_path"],
        out_path="results/qc/busco/unpolished_{assembler}",
        out_name="{sample}"
    container: config["containers"]["busco"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_path}"
        busco -i "{input.assembly}" \
              -l "{params.lineage}" \
              -o "{params.out_name}" \
              --out_path "{params.out_path}" \
              --download_path "{params.db_path}" \
              -m genome \
              --cpu {threads} \
              --force
        mkdir -p "$(dirname {output.summary})"
        find "{params.out_path}/{params.out_name}" -name "short_summary*.txt" -exec cp {{}} "{output.summary}" \;
        test -s "{output.summary}"
        """

# quast_unpolished_standard: QUAST contiguity stats, contigs >= 1 kb
rule quast_unpolished_standard:
    input:
        assembly="results/assembly/{assembler}/{sample}/assembly.fasta"
    output:
        report_txt="results/qc/quast/unpolished_{assembler}/{sample}/report.txt",
        report_html="results/qc/quast/unpolished_{assembler}/{sample}/report.html"
    params:
        outdir="results/qc/quast/unpolished_{assembler}/{sample}",
        min_contig=config["qc"]["quast"]["min_contig_standard"]
    container: config["containers"]["quast"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        quast.py "{input.assembly}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --min-contig {params.min_contig} \
            --fungus
        test -s "{output.report_txt}"
        """

# ============================================================================
# DIAGNOSTIC QC: mapping + barrnap + quast
# ============================================================================
# map_reads_to_assembly: minimap2 read-back alignment for QC and downstream tracks
rule map_reads_to_assembly:
    input:
        assembly=ref_fasta,
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        bam="results/qc/mapping/{sample}/{sample}.MD.sorted.bam",
        bai="results/qc/mapping/{sample}/{sample}.MD.sorted.bam.bai"
    container: config["containers"]["medaka"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})"
        minimap2 -t {threads} -ax map-ont --MD "{input.assembly}" "{input.reads}" | \
            samtools sort -@ {threads} -o "{output.bam}"
        samtools index -@ {threads} "{output.bam}"
        test -s "{output.bam}"
        test -s "{output.bai}"
        """

# assembly_mapping_stats: flagstat, per-contig coverage, and samtools stats
rule assembly_mapping_stats:
    input:
        bam="results/qc/mapping/{sample}/{sample}.MD.sorted.bam",
        bai="results/qc/mapping/{sample}/{sample}.MD.sorted.bam.bai"
    output:
        flagstat="results/qc/mapping/{sample}/{sample}.flagstat.txt",
        coverage="results/qc/mapping/{sample}/{sample}.coverage.tsv",
        stats="results/qc/mapping/{sample}/{sample}.stats.txt"
    container: config["containers"]["medaka"] 
    threads: 4
    shell:
        r"""
        set -euo pipefail
        samtools flagstat -@ {threads} "{input.bam}" > "{output.flagstat}"
        samtools coverage -H "{input.bam}" > "{output.coverage}"
        samtools stats -@ {threads} "{input.bam}" > "{output.stats}"
        test -s "{output.stats}"
        """

# barrnap_rrna: rRNA gene annotation; used to flag rDNA-associated regions
rule barrnap_rrna:
    input:
        assembly=ref_fasta,
    output:
        gff="results/qc/barrnap/{sample}/barrnap.gff",
        contigs="results/qc/barrnap/{sample}/rRNA_contigs.txt"
    container: config["containers"]["barrnap"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gff})"
        barrnap --kingdom euk --threads {threads} "{input.assembly}" > "{output.gff}"
        awk '$3=="rRNA" {{print $1}}' "{output.gff}" | sort -u > "{output.contigs}"
        touch "{output.contigs}"
        """

# quast_core: QUAST restricted to large contigs (>= 1 Mb)
rule quast_core:
    input:
        assembly="results/medaka/{sample}/consensus.fasta"
    output:
        report_txt="results/qc/quast_core/{sample}/report.txt"
    params:
        outdir="results/qc/quast_core/{sample}",
        min_contig=config["qc"]["quast"]["min_contig_core"]
    container: config["containers"]["quast"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        quast.py "{input.assembly}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --min-contig {params.min_contig} \
            --fungus
        test -s "{output.report_txt}"
        """

# ============================================================================
# POLISHING
# ============================================================================
def racon_ref(wildcards):
    if int(wildcards.round) == 1:
        return f"results/assembly/hifiasm/{wildcards.sample}/assembly.fasta"
    else:
        prev_round = int(wildcards.round) - 1
        return f"results/racon/{wildcards.sample}/racon_r{prev_round}.fasta"

# racon_map: minimap2 alignment used as input for one Racon round
rule racon_map:
    input:
        ref=racon_ref,
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        paf="results/racon/{sample}/round{round}.paf"
    threads: config["threads"]["medium"]
    container: config["containers"]["minimap2"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.paf})"
        zcat "{input.reads}" | minimap2 -t {threads} -x map-ont "{input.ref}" - > "{output.paf}"
        test -s "{output.paf}"
        """

# racon_polish: one round of Racon consensus correction
rule racon_polish:
    input:
        ref=racon_ref,
        reads="results/filtered/{sample}_clean.fastq.gz",
        paf="results/racon/{sample}/round{round}.paf"
    output:
        fasta="results/racon/{sample}/racon_r{round}.fasta"
    threads: config["threads"]["medium"]
    container: config["containers"]["racon"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta})"
        racon -t {threads} "{input.reads}" "{input.paf}" "{input.ref}" > "{output.fasta}"
        test -s "{output.fasta}"
        """

# medaka: final polishing with Medaka, model matched to basecalling chemistry
rule medaka:
    input:
        assembly=lambda w: f"results/racon/{w.sample}/racon_r{config['polishing']['racon_rounds']}.fasta",
        reads="results/filtered/{sample}_clean.fastq.gz"
    output:
        fasta="results/medaka/{sample}/consensus.fasta"
    params:
        outdir="results/medaka/{sample}",
        model=config["polishing"]["medaka_model"]
    threads: config["threads"]["medium"]
    container: config["containers"]["medaka"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        tmpfq=$(mktemp)
        zcat "{input.reads}" > "$tmpfq"
        medaka_consensus -i "$tmpfq" -d "{input.assembly}" -o "{params.outdir}" -m "{params.model}" -t {threads}
        rm -f "$tmpfq"
        test -s "{output.fasta}"
        """

# ============================================================================
# POLISHING QC (QUAST + BUSCO)
# ============================================================================
# quast_racon: QUAST after each Racon round, used to track polishing effects
rule quast_racon:
    input:
        assembly="results/racon/{sample}/racon_r{round}.fasta"
    output:
        report_txt="results/qc/quast_racon/{sample}_r{round}/report.txt"
    params:
        outdir="results/qc/quast_racon/{sample}_r{round}",
        min_contig=config["qc"]["quast"]["min_contig_standard"]
    container: config["containers"]["quast"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        quast.py "{input.assembly}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --min-contig {params.min_contig} \
            --fungus
        test -s "{output.report_txt}"
        """

# quast_polished: QUAST after final Medaka polishing
rule quast_polished:
    input:
        assembly="results/medaka/{sample}/consensus.fasta"
    output:
        report_txt="results/qc/quast_polished/{sample}/report.txt"
    params:
        outdir="results/qc/quast_polished/{sample}",
        min_contig=config["qc"]["quast"]["min_contig_standard"]
    container: config["containers"]["quast"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        quast.py "{input.assembly}" \
            -o "{params.outdir}" \
            --threads {threads} \
            --min-contig {params.min_contig} \
            --fungus
        test -s "{output.report_txt}"
        """

# busco_racon: BUSCO after each Racon round
rule busco_racon:
    input:
        assembly="results/racon/{sample}/racon_r{round}.fasta"
    output:
        summary="results/qc/busco_racon/{sample}/r{round}/short_summary.txt"
    params:
        lineage=config["qc"]["busco_lineage"],
        db_path=config["qc"]["busco_db_path"],
        out_path="results/qc/busco_racon/{sample}/r{round}",
        out_name=lambda w: f"{w.sample}_r{w.round}"
    container: config["containers"]["busco"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_path}"
        busco -i "{input.assembly}" \
              -l "{params.lineage}" \
              -o "{params.out_name}" \
              --out_path "{params.out_path}" \
              --download_path "{params.db_path}" \
              -m genome \
              --cpu {threads} \
              --force
        mkdir -p "$(dirname {output.summary})"
        find "{params.out_path}/{params.out_name}" -name "short_summary*.txt" -exec cp {{}} "{output.summary}" \;
        test -s "{output.summary}"
        """

# busco_polished: BUSCO on the final polished assembly
rule busco_polished:
    input:
        assembly="results/medaka/{sample}/consensus.fasta"
    output:
        summary="results/qc/busco_polished/{sample}/short_summary.txt"
    params:
        lineage=config["qc"]["busco_lineage"],
        db_path=config["qc"]["busco_db_path"],
        out_path="results/qc/busco_polished/{sample}",
        out_name=lambda w: f"{w.sample}_medaka"
    container: config["containers"]["busco"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_path}"
        busco -i "{input.assembly}" \
              -l "{params.lineage}" \
              -o "{params.out_name}" \
              --out_path "{params.out_path}" \
              --download_path "{params.db_path}" \
              -m genome \
              --cpu {threads} \
              --force
        mkdir -p "$(dirname {output.summary})"
        find "{params.out_path}/{params.out_name}" -name "short_summary*.txt" -exec cp {{}} "{output.summary}" \;
        test -s "{output.summary}"
        """

# ============================================================================
# HETEROZYGOSITY EXPLORATION (on consensus assembly)
# ============================================================================
# clair3_align_hifiasm: align variant-calling reads back to the consensus assembly
rule clair3_align_hifiasm:
    input:
        assembly=ref_fasta,
        reads="results/filtered_phasing/{sample}_clean.fastq.gz"
    output:
        bam="results/clair3_hifiasm/{sample}/reads_to_hifiasm.bam",
        bai="results/clair3_hifiasm/{sample}/reads_to_hifiasm.bam.bai"
    container:
        config["containers"]["medaka"]
    threads:
        config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})

        minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} \
          | samtools sort -@ {threads} -o {output.bam}

        samtools index -@ {threads} {output.bam}

        test -s {output.bam}
        test -s {output.bai}
        """

# clair3_call_hifiasm: variant calling with Clair3 (ONT R10.4.1 super-accurate model)
rule clair3_call_hifiasm:
    input:
        assembly=ref_fasta,
        bam="results/clair3_hifiasm/{sample}/reads_to_hifiasm.bam",
        bai="results/clair3_hifiasm/{sample}/reads_to_hifiasm.bam.bai"
    output:
        vcf="results/clair3_hifiasm/{sample}/merge_output.vcf.gz",
        tbi="results/clair3_hifiasm/{sample}/merge_output.vcf.gz.tbi"
    params:
        outdir="results/clair3_hifiasm/{sample}",
        model_path=config["clair3_model"]
    container: config["containers"]["clair3"]
    threads: config["threads"]["medium"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        if [ ! -f "{input.assembly}.fai" ]; then samtools faidx "{input.assembly}"; fi
        run_clair3.sh --bam_fn="{input.bam}" --ref_fn="{input.assembly}" --output="{params.outdir}" \
            --threads={threads} --platform="ont" --model_path="{params.model_path}" \
            --sample_name="{wildcards.sample}" --include_all_ctgs
        test -s "{output.vcf}"
        """

# clair3_het_stats_hifiasm: genome-wide heterozygosity statistics from PASS variants
rule clair3_het_stats_hifiasm:
    input:
        vcf="results/clair3_hifiasm/{sample}/merge_output.vcf.gz",
        assembly=ref_fasta
    output:
        stats="results/clair3_hifiasm/{sample}/het_stats.txt"
    container: config["containers"]["bcftools"]
    threads: 1
    shell:
        r"""
        set -euo pipefail

        # index fasta if needed
        if [ ! -f "{input.assembly}.fai" ]; then
            samtools faidx "{input.assembly}"
        fi

        ASM_BP=$(awk '{{s+=$2}} END {{print s}}' "{input.assembly}.fai")

        TOTAL=$(bcftools view -H -f PASS "{input.vcf}" | wc -l || true)
        HET=$(bcftools view -f PASS "{input.vcf}" | bcftools query -f '[%GT]\n' | grep -c "0/1" || true)
        HOM=$(bcftools view -f PASS "{input.vcf}" | bcftools query -f '[%GT]\n' | grep -c "1/1" || true)

        {{
          echo "=== Heterozygosity stats for {wildcards.sample} ==="
          echo "Total PASS variants: $TOTAL"
          echo "Het variants (0/1): $HET"
          echo "Hom-alt variants (1/1): $HOM"
          echo "Assembly size (bp): $ASM_BP"
          if [ "$HET" -gt 0 ]; then
            awk -v het="$HET" -v asm="$ASM_BP" 'BEGIN {{printf "Heterozygous variants per Mb: %.2f\n", het/(asm/1000000)}}'
          fi
        }} > "{output.stats}"

        test -s "{output.stats}"
        """

# het_density_clair3: per-contig and 100 kb windowed heterozygosity density
rule het_density_clair3:
    input:
        vcf="results/clair3_hifiasm/{sample}/merge_output.vcf.gz",
        assembly=ref_fasta
    output:
        contig_tsv="results/qc/het_density_hifiasm/{sample}/contig_het_density.tsv",
        windows_tsv="results/qc/het_density_hifiasm/{sample}/het_windows.tsv",
        plot_contig="results/qc/het_density_hifiasm/{sample}/het_density_by_contig.png",
        plot_windows="results/qc/het_density_hifiasm/{sample}/het_density_windows_topcontigs.png"
    params:
        window_size=100000,
        top_n=12
    conda:
        "envs/hetplot.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.contig_tsv})"
        if [ ! -f "{input.assembly}.fai" ]; then
            samtools faidx "{input.assembly}"
        fi
        export MPLBACKEND=Agg
        python3 scripts/het_density.py \
            --vcf "{input.vcf}" \
            --fai "{input.assembly}.fai" \
            --window {params.window_size} \
            --top {params.top_n} \
            --out_contigs "{output.contig_tsv}" \
            --out_windows "{output.windows_tsv}" \
            --plot_contigs "{output.plot_contig}" \
            --plot_windows "{output.plot_windows}"
        """

# ============================================================================
# TELOMERE MOTIF EXPLORATION & CIRCOS PLOT INPUT DATA GENERATION
# ============================================================================
# genome_windows: 10 kb non-overlapping windows over the assembly
rule genome_windows:
    input:
        fasta="results/medaka/{sample}/consensus.fasta"
    output:
        windows = f"results/circos/{{sample}}/genome_windows_10kb.bed"
    container: config["containers"]["bedtools"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.windows})

        # Create .fai if needed
        if [ ! -f "{input.fasta}.fai" ]; then
            samtools faidx "{input.fasta}"
        fi

        cut -f1,2 "{input.fasta}.fai" > {output.windows}.genome
        bedtools makewindows -g {output.windows}.genome -w 10000 \
            > {output.windows}
        rm -f {output.windows}.genome
        test -s {output.windows}
        """

# mosdepth_coverage: per-window read depth in 10 kb bins
rule mosdepth_coverage:
    input:
        bam="results/qc/mapping/{sample}/{sample}.MD.sorted.bam",
        bai="results/qc/mapping/{sample}/{sample}.MD.sorted.bam.bai"
    output:
        bed = f"results/circos/{{sample}}/coverage.regions.bed.gz"
    params:
        prefix = f"results/circos/{{sample}}/coverage"
    container: config["containers"]["mosdepth"]
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})

        mosdepth --by 10000 --no-per-base --threads {threads} \
            {params.prefix} {input.bam}

        test -s {output.bed}
        """

# gc_content: GC fraction per 10 kb window
rule gc_content:
    input:
        fasta="results/medaka/{sample}/consensus.fasta",
        windows = f"results/circos/{{sample}}/genome_windows_10kb.bed"
    output:
        gc = f"results/circos/{{sample}}/gc_content_10kb.tsv"
    container: config["containers"]["bedtools"]
    shell:
        r"""
        set -euo pipefail

        bedtools nuc -fi {input.fasta} -bed {input.windows} \
            | awk 'BEGIN{{OFS="\t"}} NR>1 {{print $1, $2, $3, $5}}' \
            > {output.gc}

        test -s {output.gc}
        """

# tidk_explore: scan for candidate telomeric repeat motifs across the assembly
rule tidk_explore:
    input:
        fasta = "results/medaka/{sample}/consensus.fasta"
    output:
        tsv = "results/circos/{sample}/tidk_explore_top.tsv"
    container: config["containers"]["tidk"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        tidk explore --minimum 5 --maximum 12 {input.fasta} \
            > {output.tsv}
        test -s {output.tsv}
        """

# tidk_search_ttaggg: search for the canonical fungal telomere TTAGGG (10 kb windows)
rule tidk_search_ttaggg:
    input:
        fasta = "results/medaka/{sample}/consensus.fasta"
    output:
        tsv = "results/circos/{sample}/tidk_ttaggg_telomeric_repeat_windows.tsv"
    params:
        motif  = "TTAGGG",
        outdir = "results/circos/{sample}",
        prefix = "tidk_ttaggg"
    container: config["containers"]["tidk"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        tidk search \
            --string {params.motif} \
            --output {params.prefix} \
            --dir {params.outdir} \
            --window 10000 \
            {input.fasta}
        test -s {output.tsv}
        """

# tidk_search_ttagg: search for the shorter telomere variant TTAGG
rule tidk_search_ttagg:
    input:
        fasta = "results/medaka/{sample}/consensus.fasta"
    output:
        tsv = "results/circos/{sample}/tidk_ttagg_telomeric_repeat_windows.tsv"
    params:
        motif  = "TTAGG",
        outdir = "results/circos/{sample}",
        prefix = "tidk_ttagg"
    container: config["containers"]["tidk"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        tidk search \
            --string {params.motif} \
            --output {params.prefix} \
            --dir {params.outdir} \
            --window 10000 \
            {input.fasta}
        test -s {output.tsv}
        """

# meryl_count_reads: k-mer counting from filtered reads (k=18)
rule meryl_count_reads:
    input:
        reads = "results/filtered/{sample}_clean.fastq.gz"
    output:
        meryl_db = directory("results/merqury/{sample}/reads.meryl")
    params:
        k = 18
    threads: 4
    container: config["containers"]["merqury"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.meryl_db})
        meryl k={params.k} threads={threads} memory=24 \
            count output {output.meryl_db} {input.reads}
        """

# merqury_eval: Merqury QV, k-mer completeness, and copy-number spectrum
rule merqury_eval:
    input:
        meryl_db = "results/merqury/{sample}/reads.meryl",
        asm      = "results/medaka/{sample}/consensus.fasta"
    output:
        qv = "results/merqury/{sample}/merqury_out.qv"
    params:
        outdir = "results/merqury/{sample}",
        prefix = "merqury_out"
    threads: 4
    container: config["containers"]["merqury"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        export MERQURY=/usr/local/share/merqury

        MERYL_ABS=$(readlink -f {input.meryl_db})
        ASM_ABS=$(readlink -f {input.asm})

        cd {params.outdir}
        merqury.sh "$MERYL_ABS" "$ASM_ABS" {params.prefix}

        test -s {params.prefix}.qv
        """
