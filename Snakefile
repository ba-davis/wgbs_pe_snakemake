
# Snakefile to analyze WGBS PE data

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz")

localrules: collect_fqc_metrics, collect_trimgalore_metrics, collect_bismark_metrics, collect_bismark_dedup_metrics, join_metrics

rule all:
    input:
        #expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        #expand("data/trimming/{sample}_val_{dir}.fq.gz", sample = SAMPLES, dir = ["1", "2"]),
        #expand("data/fastqc/trim/{sample}_val_{dir}_fastqc.zip", sample = SAMPLES, dir = ["1", "2"]),
        #expand("data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam", sample = SAMPLES),
        #"data/fastqc/raw/fqc_stats.table.txt",
        #"data/trimming/trimgalore_stats.txt",
        #"data/bismark_aln/bismark_stats.txt",
        #expand("data/bismark_aln/dedup/{sample}_val_1_bismark_bt2_pe.deduplicated.bam", sample = SAMPLES),
        #"data/bismark_aln/dedup/bismark_dedup_stats.txt",
	#"data/preprocessing_metrics/metrics.txt",
        expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz", sample = SAMPLES),
	"data/ide/ide_complete.txt"
        


rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule trim_galore:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        "data/trimming/{sample}_R2.fastq.gz_trimming_report.txt",
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    conda:
        "envs/trimgalore.yaml"
    params:
        basename = "{sample}",
        outdir = "data/trimming"
    shell:
        "trim_galore --paired --basename {params.basename} -o {params.outdir} {input.fwd} {input.rev}"

rule fastqc_trim:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        fwd = "data/fastqc/trim/{sample}_val_1_fastqc.zip",
        rev = "data/fastqc/trim/{sample}_val_2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/trim"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule bismark_aln:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        "data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt",
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
	    outdir = "data/bismark_aln"
    shell:
        "bismark -p 4 {params.genome_dir} -1 {input.fwd} -2 {input.rev} -o {params.outdir} --bam"

rule collect_fqc_metrics:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/raw/fqc_stats.table.txt"
    params:
        inpath = "data/fastqc/raw"
    shell:
        "scripts/collect_fastqc_metrics_PE.sh {params.inpath}"

rule collect_trimgalore_metrics:
    input:
        expand("data/trimming/{sample}_R2.fastq.gz_trimming_report.txt", sample = SAMPLES)
    output:
        "data/trimming/trimgalore_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/trimming",
        outfile = "data/trimming/trimgalore_stats.txt"
    shell:
        "python scripts/parse.trimgalore.wgbs.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule collect_bismark_metrics:
    input:
        expand("data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt", sample = SAMPLES)
    output:
        "data/bismark_aln/bismark_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/bismark_aln",
        outfile = "data/bismark_aln/bismark_stats.txt"
    shell:
        "python scripts/parse.bismark.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule bismark_dedup:
    input:
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    output:
        "data/bismark_aln/dedup/{sample}_val_1_bismark_bt2_pe.deduplication_report.txt",
        dedup_bam = "data/bismark_aln/dedup/{sample}_val_1_bismark_bt2_pe.deduplicated.bam"
    conda:
        "envs/bismark.yaml"
    params:
        outdir = "data/bismark_aln/dedup/"
    shell:
        "deduplicate_bismark -p --output_dir {params.outdir} --bam {input.bam}"

rule collect_bismark_dedup_metrics:
    input:
        expand("data/bismark_aln/dedup/{sample}_val_1_bismark_bt2_pe.deduplication_report.txt", sample = SAMPLES)
    output:
        "data/bismark_aln/dedup/bismark_dedup_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/bismark_aln/dedup",
        outfile = "data/bismark_aln/dedup/bismark_dedup_stats.txt"
    shell:
        "python scripts/parse.bismark_dedup.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule join_metrics:
    input:
        fqc = "data/fastqc/raw/fqc_stats.table.txt",
	trim = "data/trimming/trimgalore_stats.txt",
	aln = "data/bismark_aln/bismark_stats.txt",
	dedup = "data/bismark_aln/dedup/bismark_dedup_stats.txt"
    output:
        "data/preprocessing_metrics/metrics.txt"
    params:
        outfile = "data/preprocessing_metrics/metrics.txt"
    shell:
        "join -t $'\t' {input.fqc} {input.trim} | join -t $'\t' - {input.aln} | join -t $'\t' - {input.dedup} > {params.outfile}"

rule meth_extract:
    input:
        dedup_bam = "data/bismark_aln/dedup/{sample}_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.M-bias.txt"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
        outdir = "data/meth_extract"
    shell:
        "bismark_methylation_extractor -p --comprehensive --ignore 2 --ignore_r2 2 --ignore_3prime_r2 2 --merge_non_CpG --bedGraph --cytosine_report --gzip --genome_folder {params.genome_dir} -o {params.outdir} {input.dedup_bam}" 

rule methylkit_ide:
    input:
        expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample = SAMPLES)
    output:
        "data/ide/ide_complete.txt"
    conda:
        "envs/methylkit.yaml"
    params:
        inpath = "data/meth_extract",
        metadata = config["metadata_file"],
        group_var = config["group_var"],
        min_cov = config["min_cov"],
        outdir = "data/ide",
        min_cov_merge = config["min_cov_merge"],
        perc_merge = config["perc_merge"],
        merge_regional = config["merge_regional"],
        min_cpg_region = config["min_cpg_region"],
        mpg = config["mpg"],
	merge_dirname = config["merge_dirname"],
	repeat_initial_ide = config["repeat_initial_ide"]
        
    shell:
        "Rscript scripts/run_methylkit_ide.R {params.inpath} {params.metadata} {params.group_var} {params.min_cov} {params.outdir} {params.min_cov_merge} {params.perc_merge} {params.merge_regional} {params.min_cpg_region} {params.mpg} {params.merge_dirname} {params.repeat_initial_ide}"
