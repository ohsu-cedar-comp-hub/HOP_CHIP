import os  
import datetime
import logging 

# 2025-08-04
# I'm changing the way we do annotations
# using table_annovar.pl instead b/c much easier to do annotations from multiple databases
# will need to change add_annotations step too 

bam_folder = config["bams"]
result_folder = config["results"]
annov_folder = config["annovar"]
master_folder = config['master_dir']
kept_rules = config["kept_rules"]


FILES = glob_wildcards(os.path.join(bam_folder, "{lib}.sort.bam"))
NAMES = FILES.lib 

print('BAM files detected : ', NAMES)
print('# of BAM files:', len(NAMES))

rule all: 
    input: 
        read_counts = expand(os.path.join(result_folder, "{lib}.readcount"), lib = NAMES), 
        w_freq = expand(os.path.join(result_folder, "{lib}.all_freq"), lib = NAMES), 
        filt = expand(os.path.join(result_folder, "{lib}_" + str(config["min_depth"]) + "min_d_" + str(config['min_vaf'] * 100) + "percentmin_vaf.freq") , lib = NAMES), 
        for_annov = expand(os.path.join(result_folder, "{lib}.freq.annovar_file"), lib = NAMES),
        table = expand(os.path.join(result_folder, "{lib}." + str(config["build"]) + "_multianno.txt"), lib = NAMES), 
        final_annotate = expand(os.path.join(result_folder, "{lib}.annotated_readcounts"), lib = NAMES), 
        full_merge = os.path.join(result_folder, "merged_final.csv"), 
        filt_merged = os.path.join(result_folder, datetime.date.today().strftime('%y%m%d') + "_merged_filt_nonsyn_" + str(config["min_vaf"] * 100) + "to" + str(config["max_vaf"] * 100) + "percent_" + str(config["min_supp_reads"]) + "supp_BQ" + str(config["minBQ"])), 
        metrics = expand(os.path.join(result_folder, "{lib}.metrics"), lib = NAMES)



rule readcount: 
    input: 
        os.path.join(bam_folder, "{lib}.sort.bam")
    output: 
        read_counts = os.path.join(result_folder, "{lib}.readcount"), 
        metrics = os.path.join(result_folder, "{lib}.metrics")
    shell: 
        """
        RULE_NAME={rule} \
        {config[readcount]} \
            --max_depth={config[max_depth]} \
            --regions={config[regions]} \
            --ref={config[ref]} \
            --in={input} \
            --metrics={output.metrics} \
            --out={output.read_counts}
        """


rule calc_freq: 
    input: 
        read_counts = os.path.join(result_folder, "{lib}.readcount")
    params:
        metrics = os.path.join(result_folder, "{lib}.metrics")
    output: 
        w_freq = os.path.join(result_folder, "{lib}.all_freq")
    shell: 
        """
        RULE_NAME={rule} \
        {config[calc_freq]} \
            --script={config[readcount2freq]} \
            --in={input.read_counts} \
            --min_depth=20 \
            --min_vaf=0.0001 \
            --supp_reads=3 \
            --metrics={params.metrics} \
            --out={output.w_freq}
        """


rule filt: 
    input: 
        w_freq = os.path.join(result_folder, "{lib}.all_freq")
    params:
        metrics = os.path.join(result_folder, "{lib}.metrics")
    output: 
        filt = os.path.join(result_folder, "{lib}_" + str(config["min_depth"]) + "min_d_" + str(config['min_vaf'] * 100) + "percentmin_vaf.freq")  
    shell: 
        """
        RULE_NAME={rule} \
        {config[filter]} \
            --in={input.w_freq} \
            --min_depth={config[min_depth]} \
            --min_vaf={config[min_vaf]} \
            --metrics={params.metrics} \
            --out={output.filt}
        """


rule format_annov: 
    input: 
       os.path.join(result_folder, "{lib}_" + str(config["min_depth"]) + "min_d_" + str(config['min_vaf'] * 100) + "percentmin_vaf.freq")   
    output: 
        for_annov = os.path.join(result_folder, "{lib}.freq.annovar_file")
    resources: 
        time ="2:00:00"
    shell: 
        """
        {config[format_annovar]} \
            --script={config[freq2bed]} \
            --in={input} \
            --out={output.for_annov}
        """

rule do_annovar: 
    input: 
        os.path.join(result_folder, "{lib}.freq.annovar_file")
    output: 
        table = os.path.join(result_folder, "{lib}." + str(config["build"]) + "_multianno.txt")
    params: 
        annovar = os.path.join(annov_folder , "table_annovar.pl") , 
        db = os.path.join(annov_folder, "humandb")
    shell: 
        """
        {config[do_annovar]} \
            --in={input} \
            --annovar={params.annovar} \
            --db={params.db} \
            --addl={config[addl_db]} \
            --build={config[build]} \
            --outdir={result_folder}
        """

rule add_annotations: 
    input: 
        freq_file = os.path.join(result_folder, "{lib}_" + str(config["min_depth"]) + "min_d_" + str(config['min_vaf'] * 100) + "percentmin_vaf.freq"), 
        table = os.path.join(result_folder, "{lib}." + str(config["build"]) + "_multianno.txt")
    params:
        metrics = os.path.join(result_folder, "{lib}.metrics")
    output: 
        final_annotate = os.path.join(result_folder,"{lib}.annotated_readcounts")
    resources: 
        time ="6:00:00"
    shell: 
        """
        RULE_NAME={rule} \
        {config[add_annotations]} \
            --in={input.freq_file} \
            --script={config[annotate_bed]} \
            --annot={input.table} \
            --metrics={params.metrics} \
            --out={output.final_annotate}
        """


#per lib, let's do a enhanced metadata file? 
# lib name, age, gener

# want to merge the annotate files when all bam files have been processed 


rule merge: 
    input: 
        annot_file = expand(os.path.join(result_folder, "{lib}.annotated_readcounts"), lib = NAMES)
    output: 
        full_merge = os.path.join(result_folder, "merged_final.csv")
    threads: 2
    shell: 
        """

        complete=true 
        for file in {input.annot_file}; do 
            if [ ! -f "$file" ]; then 
                complete=false 
                break 
            fi 
        done 


        if [ "$complete" = true ]; then 
            srun mkdir-scratch.sh
            SCRATCH_PATH="/mnt/scratch/${{SLURM_JOB_ID}}"
            for file in {input.annot_file}; do 
                cp $file $SCRATCH_PATH
            done

            cat $SCRATCH_PATH/*.annotated_readcounts | grep -v library > $SCRATCH_PATH/merged_no_header.csv

            HEADER_FILE=$(ls $SCRATCH_PATH/*.annotated_readcounts | head -n 1)
            
            head -1 $HEADER_FILE | cat - $SCRATCH_PATH/merged_no_header.csv > $SCRATCH_PATH/merged.csv 

            mv $SCRATCH_PATH/merged.csv {output.full_merge}

            srun rmdir-scratch.sh 
        
        fi 


        """

rule final_filt: 
    input: 
        all_merged = os.path.join(result_folder, "merged_final.csv") 
    output: 
        filt_merged = os.path.join(result_folder, datetime.date.today().strftime('%y%m%d') + "_merged_filt_nonsyn_" + str(config["min_vaf"] * 100) + "to" + str(config["max_vaf"] * 100) + "percent_" + str(config["min_supp_reads"]) + "supp_BQ" + str(config["minBQ"]))
    threads: 2
    shell: 
        """
        srun mkdir-scratch.sh

        SCRATCH_PATH="/mnt/scratch/${{SLURM_JOB_ID}}"

        cp {input.all_merged} $SCRATCH_PATH

        sed 's/\t\t/\t/g' $SCRATCH_PATH/$(basename {input.all_merged}) > $SCRATCH_PATH/rm_extra_tabs.csv
        
        awk -F'\t' '$11 == "nonsynonymous" && $8 >= {config[min_vaf]} && $8 <= {config[max_vaf]} && $9 >= {config[minBQ]} && $7 >= {config[min_supp_reads]}' $SCRATCH_PATH/rm_extra_tabs.csv > $SCRATCH_PATH/output.csv

        mv $SCRATCH_PATH/output.csv {output.filt_merged}

        srun rmdir-scratch.sh 

        """


onsuccess: 

    def cleanup():
        log_path = os.path.join(result_folder, 'cleaned.log_file')
        with open(log_path, 'w') as logfile:

            if not config.get("clean", False):
                logfile.write("\nSkipped cleanup because 'clean = False' in config.\n")
                return

            logfile.write("\nStarting cleanup process...\n")

            all_outputs = set()
            kept_outputs = set()

            for step in workflow.rules:
                for f in step.output:
                    all_outputs.add(str(f))

            for name in kept_rules:
                step = workflow.get_rule(name)
                for f in step.output:
                    kept_outputs.add(str(f))

            tbd_removed = all_outputs - kept_outputs
            extensions = {os.path.splitext(name)[1] for name in tbd_removed if '.' in name}

            for filename in os.listdir(result_folder):
                if any(filename.endswith(ext) for ext in extensions):
                    filepath = os.path.join(result_folder, filename)
                    if os.path.isfile(filepath):
                       logfile.write('\nRemoving file:' + str(filepath) + '\n')
                       os.remove(filepath)  

    cleanup()

   







