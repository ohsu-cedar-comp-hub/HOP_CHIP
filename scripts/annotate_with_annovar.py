# import standard libraries
import numpy as np
import pandas as pd
import sys
import os
from argparse import ArgumentParser
import json
from datetime import datetime
from collections import Counter
import time
import re
import glob
from decimal import Decimal
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main():
    #Parameters to be input.
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="", required=True)
    parser.add_argument("--outdir", action="store", dest="outdir", help="", default='./')
    parser.add_argument("--annov", action="store", dest="annov", help="")
    o = parser.parse_args()

    in_calls = o.infile #*.freq

    in_annovar = o.annov # *.txt
    outdir = o.outdir

    calls = pd.read_csv(in_calls, sep='\t') 
    annovar = pd.read_csv(in_annovar, sep='\t')
    lib = calls.iloc[0,0]
    # filter only extronic 
    exonic_only = annovar[annovar['ExonicFunc.refGene'] != '.']

    print('Annotating calls with Annovar output...')
    print('Input freq files:', in_calls)
    print('Input Annovar file:', in_annovar)
    print(f'Outfile: {outdir}/{lib}.annotated_readcounts')

    # we are going to join via key CHRPOSALT 
    exonic_only.loc[:, 'CHRPOSALT'] = exonic_only.apply(
        lambda row: f"{row.Chr}:{row.Start}:{row.Alt}" if '-' not in row.Ref 
        else f"{row.Chr}:{row.Start}:+{row.Alt}",
        axis=1
    )
    exonic_only = exonic_only.drop_duplicates(subset='CHRPOSALT', keep='first')

    calls['CHRPOSALT'] = calls.apply(lambda row: str(row.chrom) + ':' + str(row.position) + ':-' if '-' in str(row.base) else str(row.chrom) + ':' + str(row.position) + ':' + str(row.base), axis = 1) 

    # create combined df - output 
    combined = pd.merge(calls,exonic_only, on=['CHRPOSALT'],how= 'left' )
    # only want rows where annotations were found 
    combined = combined[combined['Func.refGene'].notnull()] 

    # organize cols in desired format 
    combined[['type', 'class']] = combined['ExonicFunc.refGene'].str.split(' ', n =1, expand = True)
    combined = combined.drop(['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'CHRPOSALT'], axis = 1)

    combined = combined.rename(columns={'chrom' : 'chr', 'position':'pos', 'base':'alt', 'avg_basequality': 'avg_BQ', 'Gene.refGene':'gene', 
        'AAChange.refGene': 'func'})

    new_order = ['library','chr','pos','ref','alt', 'depth', 'count', 'vaf', 'avg_BQ', 'avg_pos_as_fraction', 'type', 'class', 'gene', 'func', 'cosmic102', 'snp138']

    combined = combined[new_order]

    combined.to_csv(f'{outdir}/{lib}.annotated_readcounts', sep='\t', index=False, header=True) 





main()


