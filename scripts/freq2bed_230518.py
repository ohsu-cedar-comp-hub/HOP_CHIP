#!/usr/bin/env python

'''format bam-readcount parser_results files for Annovar (normal bed file, chr start stop ref alt)

This script is required to handle differences in indel annotation. Bam-readcount parser shows insertions as C +A (for C > CA) and Annovar expects "- A" but I'm not sure if the position needs to be moved ahead one base...

'''

import pandas as pd
import os
import sys
from argparse import ArgumentParser
from timeit import default_timer as timer
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

#from py_imports import *

def main():
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="REQUIRED, input bam-readcount parser results file ", default='sys.stdin', required=True)
    o = parser.parse_args()
    input_df = pd.read_csv(o.infile, sep='\t')
    annovar_fmt_df = pd.DataFrame([], columns=['CHR', 'start', 'stop', 'ref', 'alt'])
    annovar_fmt_dict = {}
    for index, row in input_df.iterrows():
        chr = row[1]
        pos = row[2]
        ref = row[3]
        alt = row[4]
        if "+" in alt: #insertion
            alt_ins = alt.strip('+')
            annovar_fmt_dict = {'CHR': chr, 'start': int(pos), 'stop': int(pos), 'ref': "-", 'alt': alt_ins}
            annovar_fmt_df = pd.concat([annovar_fmt_df, pd.DataFrame([annovar_fmt_dict])], ignore_index=True)
        elif "-" in alt:
            annovar_fmt_df = pd.concat([annovar_fmt_df, pd.DataFrame([annovar_fmt_dict])], ignore_index=True)
            alt_del = alt.strip('-')
            stop_adj = len(alt_del)-1            
            annovar_fmt_dict = {'CHR': chr, 'start': int(pos), 'stop': int(pos) + stop_adj, 'ref': alt_del, 'alt': '-'}
            annovar_fmt_df = pd.concat([annovar_fmt_df,pd.DataFrame([annovar_fmt_dict])], ignore_index=True)
        else:
            annovar_fmt_dict = {'CHR': chr, 'start': int(pos), 'stop': int(pos), 'ref': ref, 'alt': alt}
            annovar_fmt_df = pd.concat([annovar_fmt_df,pd.DataFrame([annovar_fmt_dict])], ignore_index=True) 
    print(annovar_fmt_df)
    annovar_fmt_df.to_csv('{}.annovar_file'.format(o.infile), mode='w', sep='\t',index=False, header=False)

main()