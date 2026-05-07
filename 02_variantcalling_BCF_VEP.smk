import os  
import datetime
import logging 
import glob


bam_folder = config["bams"]
result_folder = config["results"]

master_folder = config['master_dir']


BAM_NAMES = glob.glob(os.path.join(bam_folder, "*.sort.bam"))

# Extract sample names 
SAMPLE_NAMES = [os.path.basename(f).replace(".sort.bam", "") for f in BAM_NAMES]

# Keep only samples that have both bam and indexed files
NAMES = [
    name for name in SAMPLE_NAMES
    if os.path.exists(os.path.join(bam_folder, f"{name}.sort.bam.bai"))
]

print('BAM files detected: ', NAMES)
print('# of BAM files:', len(NAMES))

rule all: 
    input: 
        vcf = expand(os.path.join(result_folder, "bcf_calls/{lib}.raw.vcf.gz"), lib = NAMES), 
        vcf_index = expand(os.path.join(result_folder, "bcf_calls/{lib}.raw.vcf.gz.csi"), lib = NAMES), 
        filt = expand(os.path.join(result_folder, "bcf_calls/{lib}_VarQ" + 
            str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz"), lib = NAMES), 
        filt_index = expand(os.path.join(result_folder, "bcf_calls/{lib}_VarQ" + 
            str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz.csi"),lib = NAMES), 
        gnomad_filt = expand(os.path.join(result_folder, "bcf_calls/{lib}_MaxAF" +
            str(config["max_af"]) + "_VarQ" + str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz"), lib = NAMES), 
        gnomad_filt_index = expand(os.path.join(result_folder, "bcf_calls/{lib}_MaxAF" +
            str(config["max_af"]) + "_VarQ" + str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz.csi"), lib = NAMES), 
        annotated = expand(os.path.join(result_folder, "annotated_vep/{lib}.csq.vcf.gz"), lib = NAMES), 
        annotated_index = expand(os.path.join(result_folder, "annotated_vep/{lib}.csq.vcf.gz.csi"), lib = NAMES)


rule bcf_call: 
    input: 
        os.path.join(bam_folder, "{lib}.sort.bam")
    params:
        bcf_dir = os.path.join(result_folder, "bcf_calls")
    output: 
        vcf = os.path.join(result_folder, "bcf_calls/{lib}.raw.vcf.gz"), 
        vcf_index = os.path.join(result_folder, "bcf_calls/{lib}.raw.vcf.gz.csi")    
    resources: 
        time ="5:00:00"
    shell: 
        """
        {config[bcf_call]} \
            --regions={config[regions]} \
            --ref={config[ref]} \
            --in={input} \
            --min_map_q={config[min_map_q]} \
            --min_base_q={config[min_base_q]} \
            --outdir={params.bcf_dir}
        """

rule bcf_filt: 
    input: 
        calls = os.path.join(result_folder, "bcf_calls/{lib}.raw.vcf.gz")
    output: 
        filt = os.path.join(result_folder, "bcf_calls/{lib}_VarQ" +
            str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz"), 
        filt_index = os.path.join(result_folder, "bcf_calls/{lib}_VarQ" +
            str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz.csi")
    shell: 
        """
        {config[bcf_filter]} \
            --in={input.calls} \
            --min_var_q={config[min_var_q]} \
            --min_depth={config[min_depth]} \
            --min_gq={config[min_gq]} \
            --out={output.filt}
        """

rule prefilt_w_gnomad: 
    input: 
        filt = os.path.join(result_folder, "bcf_calls/{lib}_VarQ" +
            str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz")
    output: 
        gnomad_filt = os.path.join(result_folder, "bcf_calls/{lib}_MaxAF" +
            str(config["max_af"]) + "_VarQ" + str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz"), 
        gnomad_filt_index = os.path.join(result_folder, "bcf_calls/{lib}_MaxAF" + 
            str(config["max_af"]) + "_VarQ" + str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz.csi")
    resources: 
        time ="5:00:00"
    shell: 
        """
        {config[prefilt_w_gnomad]} \
            --a={config[gnomad_ref]} \
            --af={config[max_af]} \
            --in={input.filt} \
            --out={output.gnomad_filt}
        """
    

rule annotate_vep: 
    input: 
        gnomad_filt = os.path.join(result_folder, "bcf_calls/{lib}_MaxAF" +
            str(config["max_af"]) + "_VarQ" + str(config["min_var_q"]) + "_" + str(config["min_depth"]) + "X_GQ" + str(config['min_gq']) + ".vcf.gz")
    params:
        annot_dir = os.path.join(result_folder, "annotated_vep")
    output: 
        annotated = os.path.join(result_folder, "annotated_vep/{lib}.csq.vcf.gz"),
        annotated_index = os.path.join(result_folder, "annotated_vep/{lib}.csq.vcf.gz.csi")
    threads: int(config.get("vep_t", 4))
    resources: 
        time ="5:00:00"
    shell: 
        """
        {config[annotate_vep]} \
            --in={input.gnomad_filt} \
            --vep_dir={config[vep_cache]} \
            --t={threads} \
            --build={config[build]} \
            --outdir={params.annot_dir} \
            --out={output.annotated}
        """



