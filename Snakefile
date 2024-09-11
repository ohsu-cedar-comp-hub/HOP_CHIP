
import os  

bam_folder = config["bams"]
result_folder = config["results"]
annov_folder = config["annovar"]
master_folder = config['master_dir']

FILES = glob_wildcards(os.path.join(bam_folder, "{lib}.sort.bam"))
NAMES = FILES.lib 

print('BAM files detected : ', NAMES)

rule all: 
    input: 
        read_counts = expand(os.path.join(result_folder, "{lib}.readcount"), lib = NAMES), 
        add_FN = expand(os.path.join(result_folder, "{lib}.FN.readcount"), lib = NAMES), 
        w_freq = expand(os.path.join(result_folder, "{lib}.all.freq"), lib = NAMES), 
        filt_freq = expand(os.path.join(result_folder, "{lib}.freq"), lib = NAMES), 
        for_annov = expand(os.path.join(result_folder, "{lib}.freq.annovar_file"), lib = NAMES),
        exonic_func = expand(os.path.join(result_folder, "{lib}.freq.annovar_file.exonic_variant_function"), lib = NAMES), 
        log = expand(os.path.join(result_folder,"{lib}.freq.annovar_file.log"), lib = NAMES) , 
        variant_func = expand(os.path.join(result_folder,"{lib}.freq.annovar_file.variant_function"), lib = NAMES), 
        final_annotate = expand(os.path.join(result_folder, "{lib}.annotated_readcounts"), lib = NAMES), 
        full_merge = os.path.join(master_folder, "merged_final.csv"), 
        filt_merged = os.path.join(master_folder, "merged_filt_nonsyn_" + str(config["min_vaf"] * 100) + "to" + str(config["max_vaf"] * 100) + "percent_" + str(config["min_supp_reads"]) + "supp_BQ" + str(config["minBQ"]))

rule readcount: 
    input: 
        os.path.join(bam_folder, "{lib}.sort.bam")
    output: 
        read_counts = os.path.join(result_folder, "{lib}.readcount")
    shell: 
        """
        bam-readcount -d {config[max_depth]} -w1 \
        -l {config[regions]} \
        -f {config[ref]} \
        {input}  > {output.read_counts}
        """

rule incl_FN: 
    input: 
        os.path.join(result_folder, "{lib}.readcount") 
    output: 
        add_FN = os.path.join(result_folder, "{lib}.FN.readcount")
    shell: 
        """
        awk '{{OFS="\\t"; print "{wildcards.lib}", $0}}' {input} | sed 's/.readcount//' > {output.add_FN}
        """


rule calc_freq: 
    input: 
        os.path.join(result_folder, "{lib}.FN.readcount") 
    output: 
        w_freq = os.path.join(result_folder, "{lib}.all.freq")
    shell: 
        """
        python -u {config[readcount2freq]} {input} \
        --min-cov {config[min_cov]} \
        --min-vaf {config[min_vaf]} \
        > {output.w_freq}
        """


rule initial_filt: 
    input: 
        os.path.join(result_folder, "{lib}.all.freq")
    output: 
        filt_freq = os.path.join(result_folder, "{lib}.freq")
    shell: 
        """
        awk '{{OFS=\"\\t\"; if ($7 > {config[min_depth]}) print $0}}' {input} > {output.filt_freq}
        """


rule format_annov: 
    input: 
        os.path.join(result_folder, "{lib}.freq")
    output: 
        for_annov = os.path.join(result_folder, "{lib}.freq.annovar_file")
    shell: 
        """
        python -u {config[freq2bed]} --infile {input} > {output.for_annov}
        """

rule do_annovar: 
    input: 
        os.path.join(result_folder, "{lib}.freq.annovar_file")
    output: 
        exonic_func = os.path.join(result_folder, "{lib}.freq.annovar_file.exonic_variant_function"), 
        log = os.path.join(result_folder,"{lib}.freq.annovar_file.log") , 
        variant_func = os.path.join(result_folder,"{lib}.freq.annovar_file.variant_function")
    params: 
        annovar = os.path.join(annov_folder , "annotate_variation.pl") , 
        db = os.path.join(annov_folder, "humandb")
    shell: 
        """
        {params.annovar} -build {config[build]} -exonsort {input} {params.db}
        """

rule add_annotations: 
    input: 
        freq_file = os.path.join(result_folder, "{lib}.freq"), 
        exonic_func = os.path.join(result_folder, "{lib}.freq.annovar_file.exonic_variant_function")
    output: 
        final_annotate = os.path.join(result_folder,"{lib}.annotated_readcounts")
    shell: 
        """
        python -u {config[annotate_final]} --infile {input.freq_file} --annov {input.exonic_func} --outdir {config[results]}
        """



# want to merge the annotate files when all bam files have been processed 


rule merge: 
    input: 
        annot_file = expand(os.path.join(result_folder, "{lib}.annotated_readcounts"), lib = NAMES)
    output: 
        full_merge = os.path.join(master_folder, "merged_final.csv")
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
            cat {input.annot_file} | grep -v library > merged_no_header.csv
            
            head -1 {input.annot_file[0]} | cat - merged_no_header.csv > {output.full_merge}

            rm merged_no_header.csv
        
        fi 

        """

rule final_filt: 
    input: 
        all_merged = os.path.join(master_folder, "merged_final.csv") 
    output: 
        filt_merged = os.path.join(master_folder, "merged_filt_nonsyn_" + str(config["min_vaf"] * 100) + "to" + str(config["max_vaf"] * 100) + "percent_" + str(config["min_supp_reads"]) + "supp_BQ" + str(config["minBQ"]))
    shell: 
        """
        sed 's/\t\t/\t/g' {input.all_merged} > rm_extra_tabs.csv
        awk -F'\t' '$11 == "nonsynonymous" && $8 > {config[min_vaf]} && $8 < {config[max_vaf]} && $9 >= {config[minBQ]} && $7 >= {config[min_supp_reads]}' rm_extra_tabs.csv > {output.filt_merged}

        rm rm_extra_tabs.csv

        """
