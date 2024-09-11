import pandas as pd
import sys
import os
import argparse
import re
import subprocess 

def extract_variant(s): 
    results = []
    for line in s.split(','):
        match = re.search(r'p\.(.*?)(?:,|$)', line)
        if match:
            results.append(match.group(1))
    
    return results


def main(): 
    parser = argparse.ArgumentParser() 
    parser.add_argument("--infile", "-i", help = "Custom calls file in bed format") 
    parser.add_argument("--annottool", "-a", help = "Path to annovar perl script", default = "/home/groups/CEDAR/chaoe/HOP_test/tools/annovar/annotate_variation.pl")
    parser.add_argument("--db", help = "Path to annovar human db", default = "/home/groups/CEDAR/chaoe/HOP_test/tools/annovar/humandb")
    parser.add_argument("--build", help = "genome build", default = "hg19"), 
    parser.add_argument("--outdir", '-o', help = "Path to directory where reformatted calls goes")

    o = parser.parse_args() 

    cmd = f"{o.annottool} -build {o.build} -exonsort {o.infile} {o.db}"
    result = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)

    calls = pd.read_csv(o.infile + '.exonic_variant_function', sep='[\t ]', names=['bed_line', 'type', 'class','info', 'chr', 'start', 'stop', 'ref', 'alt'], engine = 'python')

    calls = calls[calls['type'] != 'unknown']

    calls.loc[calls['type'].str.contains('stop', case=False, na=False), 'info'] = calls['class']

    calls['gene'] = calls['info'].str.split(':').str[0]

    calls['variant'] = calls['info'].apply(extract_variant)

    new_format = calls[['gene', 'variant']].explode('variant').drop_duplicates().reset_index(drop=True)

    if o.outdir is None: 
        out_path = o.infile + '_reformatted.txt'
        new_format.to_csv(o.infile + '_reformatted.txt', sep='\t', index=False, header=True)
    
    else: 
        out_path = o.outdir + '/' + os.path.basename(o.infile) + '_reformatted.txt'
        new_format.to_csv(o.outdir + '/' + os.path.basename(o.infile) + '_reformatted.txt', sep='\t', index=False, header=True)


    print("Reformatted calls saved to", out_path)
main() 