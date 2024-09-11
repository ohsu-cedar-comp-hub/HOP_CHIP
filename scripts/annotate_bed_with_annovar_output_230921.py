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

    #in_calls = '/home/groups/CEDAR/boniface/HOP_data/BAM/22KD-115F0022-1.freq'
    in_calls = o.infile

    in_annovar = o.annov 
    outdir = o.outdir

    calls = pd.read_csv(in_calls, sep='\t') 
    annovar = pd.read_csv(in_annovar, sep='[\t ]', names=['bed_line', 'type', 'class','info', 'chr', 'start', 'stop', 'ref', 'alt'])
    annovar['CHRPOSALT'] = annovar.apply(lambda row: str(row.chr) + ':' + str(row.start) + ':' + str(row.alt) if '-' not in row.ref else str(row.chr) + ':' + str(row.start) + ':+' + str(row.alt), axis=1)
    annotated_calls = annovar['CHRPOSALT'].tolist()
    #print('annotated_calls',annotated_calls)
    annovar = annovar.set_index(['CHRPOSALT'])
    annovar = annovar[~annovar.index.duplicated(keep='first')]
    lib = calls.iloc[0,0]

    columns=['library', 'chr', 'pos', 'ref', 'alt','depth', 'count',  'vaf', 'avg_BQ', 'avg_pos_as_fraction', 'type', 'class', 'gene', 'func']
    annot_calls = pd.DataFrame([],columns=columns) #initialize output data frame
    #annot_calls.to_csv(f'{outdir}/{lib}.annotated_readcounts', mode='w', sep='\t', index=False, header=True)
    annot_calls.to_csv(f'{outdir}/{lib}.annotated_readcounts', mode='w', sep='\t', index=False, header=True)

    print('Annotating calls with Annovar output...')
    print('Input freq files:', in_calls)
    print('Input Annovar file:', in_annovar)
    print(f'Outfile: {outdir}/{lib}.annotated_readcounts')

    for index, call in calls.iterrows():
        CHR = call[1]
        POS = call[2]
        REF = call[3]
        ALT = call[4]
        depth = call[6]
        vaf = call[5]
        count = call[7]
        avg_BQ = call[8]
        avg_pos_as_fraction = call[9]
        CHRPOS = str(CHR) + ":" + str(POS)
        if '-' in ALT:
            CHRPOSALT = CHRPOS + ":-"
        else:
            CHRPOSALT = CHRPOS + ":" + ALT
        if int(POS) == 25457160:
            print(f'Found {POS} (CHRPOSALT:{CHRPOSALT}) matching annovar info = {annovar.loc[CHRPOSALT]}')
        if CHRPOSALT not in annotated_calls: 
            #print(f'Could not find mutation {CHRPOSALT} in provided annotation file... May be incorrectly formatted indel.')
            continue
        annovar_info = annovar.loc[CHRPOSALT]
        type_i = annovar_info[1]
        class_i = annovar_info[2]
        gene = annovar_info[3].split(':')[0]
        func =  annovar_info[3].split(' ')[0]

        dict1={'library': lib, 'chr': CHR, 'pos': POS, 'ref': REF, 'alt': ALT, 'depth': depth, 
                'count': count, 'vaf': vaf, 'avg_BQ' : avg_BQ, 'avg_pos_as_fraction': avg_pos_as_fraction, 
                'type': type_i, 'class' : class_i, 'gene' : gene, 'details': func}
        annot_calls=pd.concat([annot_calls, pd.DataFrame([dict1])], ignore_index=True)

    annot_calls.to_csv(f'{outdir}/{lib}.annotated_readcounts', mode='a', sep='\t',index=False, header=False) 

main()