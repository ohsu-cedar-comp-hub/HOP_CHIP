import os  
import datetime
import logging 
import glob


fastq_dir = config['fastq_dir']
result_dir = config['results_dir']

R1_NAMES = glob.glob(os.path.join(fastq_dir, "*_1.fq.gz"))

# Extract sample names (without _R1.fastq.gz)
SAMPLE_NAMES = [os.path.basename(f).replace("_1.fq.gz", "") for f in R1_NAMES]

# Keep only samples that have both _R1 and _R2 files
NAMES = [
    name for name in SAMPLE_NAMES
    if os.path.exists(os.path.join(fastq_dir, f"{name}_2.fq.gz"))
]

print("Sample names are:", NAMES)
print('# of Samples Detected:', len(NAMES))

rule all: 
    input: 
        r1_report = expand(os.path.join(result_dir, "FASTQC/{lib}_1_fastqc.html"), lib = NAMES),
        r2_report = expand(os.path.join(result_dir, "FASTQC/{lib}_2_fastqc.html"), lib = NAMES), 
        screen_report = expand(os.path.join(result_dir, "FASTQScreen/{lib}_1_screen.html"), lib = NAMES), 
        r1_trimmed = expand(os.path.join(result_dir, "Trimmed/{lib}_trimmed_1.fq.gz"), lib = NAMES), 
        r2_trimmed = expand(os.path.join(result_dir, "Trimmed/{lib}_trimmed_2.fq.gz"), lib = NAMES),
        aligned = expand(os.path.join(result_dir, "Aligned/{lib}.aligned.bam"), lib = NAMES), 
        rg_checked = expand(os.path.join(result_dir, "Aligned/{lib}.rg_exists.txt"), lib = NAMES), 
        deduped = expand(os.path.join(result_dir, "Deduped/{lib}.sort.bam"), lib = NAMES), 
        deduped_index = expand(os.path.join(result_dir, "Deduped/{lib}.sort.bam.bai"), lib = NAMES)



rule fastqc: 
    input: 
        r1 = os.path.join(fastq_dir, "{lib}_1.fq.gz"), 
        r2 = os.path.join(fastq_dir, "{lib}_2.fq.gz")
    output: 
        r1_report = os.path.join(result_dir, "FASTQC/{lib}_1_fastqc.html"), 
        r2_report = os.path.join(result_dir, "FASTQC/{lib}_2_fastqc.html"), 
        r1_zip = os.path.join(result_dir, "FASTQC/{lib}_1_fastqc.zip"), 
        r2_zip = os.path.join(result_dir, "FASTQC/{lib}_2_fastqc.zip") 
    params: 
        fastqc_dir = os.path.join(result_dir, "FASTQC")
    threads: 2
    resources: 
        time ="5:00:00"
    shell: 
        """
        {config[fastqc]} \
        --r1={input.r1} \
        --r2={input.r2} \
        --outdir={params.fastqc_dir} \
        --threads={threads}
        """


rule fastqscreen: 
    input: 
        r1 = os.path.join(fastq_dir, "{lib}_1.fq.gz"), 
        r2 = os.path.join(fastq_dir, "{lib}_2.fq.gz")
    output: 
        screen_report = os.path.join(result_dir, "FASTQScreen/{lib}_1_screen.html")
    params: 
        fastqscreen_dir = os.path.join(result_dir, "FASTQScreen")
    threads: 8
    resources: 
        time ="5:00:00"
    shell: 
        """
        {config[fastq_screen]} \
        --r1={input.r1} \
        --r2={input.r2} \
        --conf={config[fastq_screen_config]} \
        --outdir={params.fastqscreen_dir} \
        --threads={threads}
        """


rule trim: 
    input: 
        r1 = os.path.join(fastq_dir, "{lib}_1.fq.gz"), 
        r2 = os.path.join(fastq_dir, "{lib}_2.fq.gz")
    output:
        trimmed_r1 = os.path.join(result_dir, "Trimmed/{lib}_trimmed_1.fq.gz"), 
        trimmed_r2 = os.path.join(result_dir, "Trimmed/{lib}_trimmed_2.fq.gz")
    params: 
        trim_dir = os.path.join(result_dir, "Trimmed")
    threads: 4
    resources: 
        time ="5:00:00", 
        mem = "20G"
    shell: 
        """
        {config[trim]} \
        --r1={input.r1} \
        --r2={input.r2} \
        --r1_adaptor={config[r1_adaptor]} \
        --r2_adaptor={config[r2_adaptor]} \
        --outdir={params.trim_dir} \
        --lib={wildcards.lib} \
        --threads={threads}
        """

rule align: 
    input:
        trimmed_r1 = os.path.join(result_dir, "Trimmed/{lib}_trimmed_1.fq.gz"), 
        trimmed_r2 = os.path.join(result_dir, "Trimmed/{lib}_trimmed_2.fq.gz")
    output: 
        aligned = os.path.join(result_dir, "Aligned/{lib}.aligned.bam") 
    params:
        align_dir = os.path.join(result_dir, "Aligned")
    threads: 10
    resources: 
        time ="18:00:00", 
        mem = "30G"
    shell: 
        """
        {config[align]} \
        --r1={input.trimmed_r1} \
        --r2={input.trimmed_r2} \
        --genome={config[genome]} \
        --outdir={params.align_dir} \
        --lib={wildcards.lib} \
        --threads={threads}
        """

rule check_rg: 
    input: 
        aligned = os.path.join(result_dir, "Aligned/{lib}.aligned.bam") 
    output: 
        rg_checked = os.path.join(result_dir, "Aligned/{lib}.rg_exists.txt")
    params: 
        align_dir = os.path.join(result_dir, "Aligned")
    resources: 
        time="8:00:00"
    shell: 
        """
        {config[check_rg]} \
        --in={input.aligned} \
        --outdir={params.align_dir} \
        --lib={wildcards.lib} \
        --rg={config[add_rg]} 
        """

rule deduplicate: 
    input: 
        aligned = os.path.join(result_dir, "Aligned/{lib}.aligned.bam"), 
        rg_checked = os.path.join(result_dir, "Aligned/{lib}.rg_exists.txt")
    output: 
        deduped = os.path.join(result_dir, "Deduped/{lib}.sort.bam"), 
        deduped_index = os.path.join(result_dir, "Deduped/{lib}.sort.bam.bai"),
        metrics = os.path.join(result_dir, "Deduped/{lib}.dedup_metrics.txt")
    params: 
        dedup_dir = os.path.join(result_dir, "Deduped")
    resources: 
        time="12:00:00", 
        mem = "20G"
    shell: 
        """
        {config[dedup]} \
        --in={input.aligned} \
        --outdir={params.dedup_dir} \
        --lib={wildcards.lib} 
        """