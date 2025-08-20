# Analyzing HOP data W/ Focus on CHIP Calls 

Given sequencing data from HOP, the user can identify CHIP calls, compare with known calls, perform various annotations and create a variety of companion visualizations. 

## Getting Started 

```bash 
git clone https://github.com/ohsu-cedar-comp-hub/HOP_CHIP.git
conda env create -f HOP.yaml 
conda activate HOP
```


### Part 1: Analyzing Bams 
#### Summary
This workflow will take your sorted indexed bams and for each, will perform read counts, calculate VAFs, filter and perform gene annotations via annovar. 
After, your results will be merged into a single dataframe for convenience and for downstream analysis. 

**Want to know how I developed this workflow?**

Go here for an in-depth guide: https://ohsuitg-my.sharepoint.com/:w:/r/personal/chaoe_ohsu_edu/Documents/Guide%20to%20Snakemake%20CHIP%20Workflow%20-%20Notes.docx?d=we5b5d1e4a231405c906d9cdc572f9541&csf=1&web=1&e=rXWT7v


Don't have access? Email me: chaoe@ohsu.edu

Tools needed: 
* Annovar (https://annovar.openbioinformatics.org/en/latest/user-guide/download/)

This is a Snakemake workflow and will run through calculating readcounts, VAF calculations, filtering,  annotating and merging in a step by step process. 

#### Running Workflow
To run this workflow, you will need preliminary knowledge of how Snakemake functions. 
Refer here for a general guide: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html

1- Start off with sorted indexed bam files (*.sorted.bam and *.sorted.bam.bai).
If you have your files handy, place them in their own directory labelled appropriately. 
If you have your files in Ceph, you will need to access them using your AWS profile. You will need to specify where you want your data to be stored and what your AWS profile is when running the launch script. Refer to Step 5. 

2- Edit config/config_bams.json accordingly and change paths as needed. 
NOTE: Change addl_db accordingly if you do not have additional databases you want annovar to annotate with. If you want more information about this go here: https://annovar.openbioinformatics.org/en/latest/user-guide/download/#additional-databases


3- Edit config/config.v8+.yaml to change cluster configuration settings. 

4- Perform a dry run on just Snakemake file to ensure no issues. 

```bash 
snakemake -n --profile simple/ --configfile=config/config_bams.json
```

5- Run workflow via launch script. Only include -d and -p if you need to pull data from Ceph. Refer to Step 1. 

```bash
sbatch AnalyzeBams.sh -c config/config_bams.json -d data/ -p xyz  
```

All output files generated for each input will be located in the results folder indicated by config/config_bams.json. 

All outputs per sample wil be: 
 *.readcount, *.all_freq, *.metrics, *.#min_d_#percentmin_vaf.freq, *.freq.annovar_file, *.hg19_multianno.txt, *.annotated_readcounts

After merging: 
merged_final.csv, 
{rundate}_merged_filt_nonsyn_{minvaf}
to{maxvaf} _ {#suppreads}supp_BQ{minBQ}


NOTE: These outputs only apply if you specified clean=False in the config/config.json. 
If you had chosen to set clean=True, you can then specify the rules in kept_rules whose outputs you want to keep. This will be entirely based on what you need the workflow for. 


NOTE: I really HIGHLY recommend refering to the in-depth guide to access more information about this workflow. There is a LOT of information in here: 

https://ohsuitg-my.sharepoint.com/:w:/r/personal/chaoe_ohsu_edu/Documents/Guide%20to%20Snakemake%20CHIP%20Workflow%20-%20Notes.docx?d=we5b5d1e4a231405c906d9cdc572f9541&csf=1&web=1&e=rXWT7v

If you need access or have any questions, please email me: chaoe@ohsu.edus


### Part 2: Annotating Counts and Analyzing Calls 
#### Summary
This script can be used in 2 ways: 

1- Directly from Part 1 by using the resulting dataframe as input and then generating a new dataframe with added annotations from various sources (ex. known calls, dbsnp, metadata). 

2- Using an already annotated dataframe and a dictionary containing the lists of known calls as inputs. 

The annotated dataframe will be used to create multiple visualizations specifically for mutation calls that match the provided known calls. 

Files Involved: 
* {Today(YYMMMDD)}_merged_filt_nonsyn_{minvaf}to{maxvaf}percent_{minreads}supp_BQ{minBQ} (from Part 1)
* scripts/AnalyzeMutCalls.py 
* config/config_annot.json 

#### Running Script 
1- Run AnalyzeMutCalls.py. 


```
usage: AnalyzeMutCalls.py [-h] [--infile INFILE] [--metadata METADATA] [--dbsnp DBSNP] [--annot ANNOT] [--plt PLT] [--ver VER] [--minvaf MINVAF] [--maxvaf MAXVAF]
                          [--mincount MINCOUNT] [--config CONFIG] [--calls KEY=VALUE [KEY=VALUE ...]] [--dict DICT]

options:
  -h, --help            show this help message and exit
  --infile INFILE, -i INFILE
                        Input File produced from Snakemake
  --metadata METADATA, -m METADATA
                        Patient metadata information
  --dbsnp DBSNP, -d DBSNP
                        Known sites in DBSNP
  --annot ANNOT, -a ANNOT
                        Annotated file
  --plt PLT, -p PLT     Directory where plots go
  --ver VER, -v VER     OPTIONAL Can specify version if multiple runs
  --minvaf MINVAF       OPTIONAL Minimum variant allele frequency
  --maxvaf MAXVAF       OPTIONAL Maximum variant allele frequency
  --mincount MINCOUNT   OPTIONAL Minimum read count
  --config CONFIG, -c CONFIG
                        OPTIONAL Path to JSON configuration file
  --calls KEY=VALUE [KEY=VALUE ...]
                        Set a number of key-value pairs where name = path to calls
  --dict DICT           Path to dictionary of calls

additional information: If you are starting with a file produced from Snakemake, --infile, --metadata, --calls, --dbsnp and --plt are required. Any that you don't
fill in will be filled in by default with config file. If you have an annotated file already, --annot, --dict and --plt are required.
```
If you are using the resulting file directly from Part 1 and you want to use your own paths to annotation files, an example would be:  

```bash 
python -u AnalyzeMutCalls.py -i {Today(YYMMMDD)}_merged_filt_nonsyn_{minvaf}to{maxvaf}percent_{minreads}supp_BQ{minBQ} -m HOP_first_500_metadata.csv -d dbsnp_22KD-115F0022-1_50X_intersect.bed -p /home/groups/CEDAR/chaoe/HOP_test/figs --calls watson=/home/groups/CEDAR/chaoe/HOP_test/calls/CHIP_sites_from_Watson_etal_Blundell_lab.txt beataml=/home/groups/CEDAR/chaoe/HOP_test/calls/CHIP_sites_from_BeatAML.txt
```

If you are using a previously annotated file and a corresponding dictionary, and just want the visualizations, an example could look like:  

```bash
python -u AnalyzeMutCalls.py -a HOP_merged_v1_240909.annotated_readcounts_crossref.csv --dict HOP_merged_v1_240909.callsdict.txt -p /home/groups/CEDAR/chaoe/HOP_test/figs
```

If your file is annotated, and you are missing the dict, you can create a .txt file and manually create one that maps each call list used in the file to the full name of the calls list. 
```
{'b':'beataml', 'c':'custom'}
```
NOTE: Because there are many parameters, you could also modify config_annot.json as an alternative to specifying parameters in the command line. 

All visualizations generated will be in their respective matched calls' subdirectories in your given plot directory. 
An example of the visualizations you will see per subdirectory is in figs/beataml of this repo. 


### Part 2a: Creating Custom Calls
#### Summary
Starting with a BED3 file, this script will convert your calls to a readable list so that it can be used above in AnalyzeMutCalls.py. 

Files Involved: 
* scripts/ConvertCalls.py 

#### Running Script 
1- Input file must be a BED3 file with format: 
```
chr start   end alt
```
2- Run ConvertCalls.py 
```
usage: ConvertCalls.py [-h] [--infile INFILE] [--annottool ANNOTTOOL] [--db DB] [--build BUILD] [--outdir OUTDIR]

options:
  -h, --help            show this help message and exit
  --infile INFILE, -i INFILE
                        Custom calls file in bed format
  --annottool ANNOTTOOL, -a ANNOTTOOL
                        Path to annovar perl script
  --db DB               Path to annovar human db
  --build BUILD         genome build
  --outdir OUTDIR, -o OUTDIR
                        Path to directory where reformatted calls goes

```
A usage example would be: 

```bash 
python -u ConvertCalls.py -i testcalls.freq.annovar_file -o /home/groups/CEDAR/chaoe/HOP_test/calls --annottool /home/groups/CEDAR/chaoe/HOP_test/tools/annovar/annotate_variation.pl --db /home/groups/CEDAR/chaoe/HOP_test/tools/annovar/humandb --build hg19

```

Resulting output is a *_reformatted.txt that will be in your specified output directory. 
This file can be added to the AnalyzeMutCalls.py as part of the --calls parameter. 

## Additional Notes 
You can enter in any number of list of calls for AnalyzeMutCalls.py. But, when using multiple list of calls, choose names that start with a different letter. The created annotation dataframe saves the first letter of the name if there is a match, so this will prevent any confusion. 
Good example: custom, beataml, watson . Bad example: custom, calls, california

