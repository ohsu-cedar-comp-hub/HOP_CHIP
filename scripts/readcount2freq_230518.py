#!/usr/bin/env python

#https://github.com/genome/bam-readcount/blob/master/tutorial/scripts/parse_brc.py

#this version has been adapted to accept files with 1 additional columns ("library")
#It is intended to convert .readcount files to .freq files
#

import os
import argparse
import logging

SCRIPT_PATH = os.path.abspath(__file__)
FORMAT = '[%(asctime)s] %(levelname)s %(message)s'
l = logging.getLogger()
lh = logging.StreamHandler()
lh.setFormatter(logging.Formatter(FORMAT))
l.addHandler(lh)
l.setLevel(logging.INFO)
debug = l.debug
info = l.info
warning = l.warning
error = l.error

DESCRIPTION = '''
Parse bam-readcount output with optional filters
'''

EPILOG = '''
'''


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass


parser = argparse.ArgumentParser(description=DESCRIPTION,
    epilog=EPILOG,
    formatter_class=CustomFormatter)

parser.add_argument('input',
    help='input file (.bam-readcount) from bam_readcount_output')
parser.add_argument('--min-cov',
    action='store',
    type=int,
    help='Minimum coverage to report variant',
    default=0)
parser.add_argument('--min-vaf',
    action='store',
    type=float,
    help='Minimum VAF to report variant',
    default=0)
parser.add_argument('--max-vaf',
    action='store',
    type=float,
    help='Maximum VAF to report variant',
    default=0.5)    
parser.add_argument('--suprds',
    action='store',
    type=int,
    help='Minimum supporting reads to report variant',
    default=1)
parser.add_argument('-v',
    '--verbose',
    action='store_true',
    help='Set logging level to DEBUG')
parser.add_argument('--omitCA',
    action='store_true',
    help='filter C to A substitutions')    

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)

debug('%s begin', SCRIPT_PATH)

headers = ['library', 'chrom', 'position', 'ref', 'base', 'vaf', 'depth', 'count',
    'avg_basequality', 'avg_pos_as_fraction'
]
print('\t'.join(headers))

# Per-base/indel data fields
# IMPORTANT: this relies on Python 3.6+ to maintain insertion order
# Each field is a key with value a function to convert to the
# appropriate data type
base_fields = {
    'base': str,
    'count': int,
    'avg_mapping_quality': float,
    'avg_basequality': float,
    'avg_se_mapping_quality': float,
    'num_plus_strand': int,
    'num_minus_strand': int,
    'avg_pos_as_fraction': float,
    'avg_num_mismatches_as_fraction': float,
    'avg_sum_mismatch_qualities': float,
    'num_q2_containing_reads': int,
    'avg_distance_to_q2_start_in_q2_reads': float,
    'avg_clipped_length': float,
    'avg_distance_to_effective_3p_end': float
}

# Open the bam-readcount output file and read it line by line
# Note that the output has no header, so we consume every line
with open(args.input) as in_fh:
  for line in in_fh:
    line = line.strip()    
    fields = line.split('\t')  
    library = fields[0]
    chrom = fields[1]              
    position = int(fields[2])      
    reference_base = fields[3]    
    depth = int(fields[4])
    for base_data_string in fields[5:]:
      # We will store per-base/indel data in a dict
      base_data = {}
      # Split the base/indel data on ':'
      base_values = base_data_string.split(':')
      # Iterate over each field of the base/indel data
      for i, base_field in enumerate(base_fields.keys()):
        # Store each field of data, converting to the appropriate
        # data type
        base_data[base_field] = base_fields[base_field](base_values[i])

      # Skip zero-depth bases
      if depth == 0:
        continue
      # Skip reference bases and bases with no counts
      if base_data['base'] == reference_base or base_data['count'] == 0 or base_data['base'] == "N":
        continue
      # Calculate an allele frequency (VAF) from the base counts
      if args.omitCA and reference_base == "C" and base_data['base'] == "A":  
        continue  
      vaf = base_data['count'] / depth
      # Filter on minimum depth and VAF
      if depth >= args.min_cov and vaf >= args.min_vaf and base_data['count'] >= args.suprds: 
        # Output count and VAF data as well as avg_pos_as_fraction
        print('\t'.join(
            str(x) for x in (library, chrom, position, reference_base, base_data['base'],
            '%0.6f' % (vaf), depth, base_data['count'],
            base_data['avg_basequality'], base_data['avg_pos_as_fraction'])))

debug('%s end', (SCRIPT_PATH))