# Analyzing HOP data W/ Focus on CHIP Calls 

Given sequencing data from HOP, the user can identify CHIP calls, compare with known calls, perform various annotations and create a variety of companion visualizations. 

## Getting Started 

```bash 
git clone https://github.com/ohsu-cedar-comp-hub/HOP_CHIP.git
conda env create -f HOP.yaml --n HOP 
conda activate HOP
```


### Part 1: Analyzing Bams 
#### Summary
This workflow will take your sorted indexed bams and for each, will perform read counts, calculate VAFs, filter and perform gene annotations via annovar. 
After, your results will be merged into a single dataframe for convenience and for downstream analysis. 

Files involved: 
* Snakemake 
* config/config_bams.json  
* simple/config.v8+.yaml, 
* AnalyzeBams.sh,
* scripts/freq2bed_230518.py,
* scripts/annotate_bed_with_annovar_output_230921.py
* scripts/readcount2freq_230518.py 

#### Running Workflow
1- Place your sorted indexed bam files (*.sorted.bam and *.sorted.bam.bai) in their own directory. 

2- Edit config_bams.json accordingly and change paths as needed.

3- Edit config.v8+.yaml to change cluster configuration settings. 

4- Perform a dry run on just Snakemake file to ensure no issues. 

```bash 
snakemake -n --configfile config/config_bams.json
```

5- Run 
```bash
sbatch AnalyzeBams.sh config/config_bams.json 
```

All output files generated for each input will be located in the results folder indicated by config_bams.json. 

Resulting final dataframe will be located in the master directory with the name following the format: 


"{Today(YYMMMDD)}_merged_filt_nonsyn_{minvaf}to{maxvaf}percent_{minreads}supp_BQ{minBQ}"

### Part 2: Annotating Counts and Analyzing Calls 
#### Summary
This script will take as input the resulting dataframe from Part 1 and generate a new dataframe with added annotations from various sources (ex. known calls, dbsnp, metadata). 

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
For example, if you are using the final file directly from Part 1 and you want to use your own paths to annotation files: 

```bash 
python -u AnalyzeMutCalls.py -i {Today(YYMMMDD)}_merged_filt_nonsyn_{minvaf}to{maxvaf}percent_{minreads}supp_BQ{minBQ} -m HOP_first_500_metadata.csv -d dbsnp_22KD-115F0022-1_50X_intersect.bed -p /home/groups/CEDAR/chaoe/HOP_test/figs --calls watson=/home/groups/CEDAR/chaoe/HOP_test/calls/CHIP_sites_from_Watson_etal_Blundell_lab.txt beataml=/home/groups/CEDAR/chaoe/HOP_test/calls/CHIP_sites_from_BeatAML.txt
```

If you are using a previously annotated file, and just want the visualizations, an example could look like:  

```bash
python -u AnalyzeMutCalls.py -a HOP_merged_v1_240909.annotated_readcounts_crossref.csv --dict HOP_merged_v1_240909.callsdict.txt -p /home/groups/CEDAR/chaoe/HOP_test/figs
```

If your file is annotated, and you are missing the dict, you can create a .txt file and manually create one that maps each call list used in the file to the full name of the calls list. 
```
{'b':'beataml', 'c':'custom'}
```

 Resulting output will be your annotated file (if you started with the file produced from Part 1), and all visualizations generated will be in their respective matched calls subdirectories in your given plot directory. 
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

## Additional Notes 
When using multiple list of calls, choose names that start with a different letter. The created annotation dataframe saves the first letter of the name if there is a match, so this will prevent any confusion. 
A good example: custom, beataml, watson . Not good : custom, calls, california

