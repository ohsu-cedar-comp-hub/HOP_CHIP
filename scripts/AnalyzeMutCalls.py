
# import standard libraries
import numpy as np
import scipy
import scipy.stats
from scipy.stats import pearsonr
import pandas as pd
import sys
import os
import argparse
import json
from datetime import datetime
#import pysam
import copy
from collections import Counter
import time
import pickle
from math import isnan
import re
import glob
from decimal import Decimal
from functools import partial
from natsort import natsorted
import textwrap 

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D 
import seaborn as sns


# Needed Functions 
def unique(list):
    unique_list = np.unique(np.array(list)) 
    return unique_list

def sort_A_using_B(A, B):
    sorted_A = [obj for _, obj in sorted(zip(B, A))].copy()
    return sorted_A

def rand_jitter(arr):
    stdev = .01 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def jitter(obj, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs):
    return plt.scatter(rand_jitter(obj), y, s=s, c=c, marker=marker, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, alpha=alpha, linewidths=linewidths, **kwargs)

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)


def validate_args(args): # validate enough args are given for script to run 
    need_annot = False 
    
    if args.infile: 
        if not all([args.metadata, args.dbsnp, args.plt, args.calls]):
            print("Error: If `--infile` is specified, all other parameters: `--metadata`, `--calls`, `--dbsnp` and `--plt` must be specified.")
            sys.exit(1)
        need_annot = True

    if args.annot:
        if not args.plt and not args.dict:
            print("Error: If `--annot` is specified, `--plt` must be specified.")
            sys.exit(1)
    
    return need_annot

def reformat(x): # reformat calls to readable form 
    calls = pd.read_csv(x, sep = '\t')
    calls['aa_pos'] = calls.apply(lambda row: re.findall(r"[^\W\d_]+|\d+", row.variant)[0] + re.findall(r"[^\W\d_]+|\d+", row.variant)[1], axis = 1)
    calls['ref_aa'] = calls.apply(lambda row: re.findall(r"[^\W\d_]+|\d+", row.variant)[0], axis = 1)
    calls['pos'] = calls.apply(lambda row: re.findall(r"[^\W\d_]+|\d+", row.variant)[1], axis = 1)
    calls['all_info'] = calls.apply(lambda row: re.findall(r"[^\W\d_]+|\d+", row.variant), axis = 1)

    calls['gene_ref_aa_pos'] =  calls.apply(lambda row: row.gene + "_" + row.aa_pos, axis = 1)

    return calls['gene_ref_aa_pos'].to_list()
   
def get_mdnmean_data(mycalls, calls, min_vaf, max_vaf, min_ct): 
    # apply additional filtering if needed 
    if min_vaf:
       calls = calls[calls['vaf'] > min_vaf] 
    if max_vaf: 
       calls = calls[calls['vaf'] < max_vaf] 
    if min_ct: 
        calls = calls[calls['count'] > min_ct]

    # Obtain unfiltered and filtered means, medians, depths, ages, vafs grouped by library 
    
    library_unfil_mdns = mycalls.groupby('library').median(numeric_only=True).sort_values('library').reset_index()
    library_unfil_means = mycalls.groupby('library').mean(numeric_only=True).sort_values('library').reset_index()
    library_mdns = calls.groupby('library').median(numeric_only=True).sort_values('library').reset_index()
    library_means = calls.groupby('library').mean(numeric_only=True).sort_values('library').reset_index()
    mdn_depths_unfil = library_unfil_mdns.depth.to_list()
    mdn_depths = library_mdns.depth.to_list()
    mdn_vafs = library_mdns.vaf.to_list()
    mdn_ages = library_mdns.age.to_list() ## mdn age is just age of the patient, but I need this to match the mdn_vafs dataset 
    mean_depth_unfil = library_unfil_means.depth.to_list()
    mean_depth = library_means.depth.to_list()
    mean_vaf = library_means.vaf.to_list()

    return calls, {
        'unfil_depths': mdn_depths_unfil,
        'depths': mdn_depths,
        'vafs': mdn_vafs,
        'ages': mdn_ages
    }, {
        'unfil_depths': mean_depth_unfil,
        'depths': mean_depth,
        'vafs': mean_vaf
    }


def get_matching_data(mycalls, calls, mds):
    lib_data = {}
    libs = []
    call_counts = []
    vafs = []
    depths = []
    ages = []
    ages_pts_having_calls = [] ## a list of ages to match the length of vafs list of lists (i.e., only ages of patients in df_hopcalls_beataml)
    gends = []
    gends_pts_having_calls = [] ## a list of depths to match the length of vafs list of lists (i.e., only depths of calls in df_hopcalls_beataml)
    genes = []

    for library in sorted(set(mycalls.library)): ## NOTE: using unfiltered library here to get call count = 0 for patients with no calls after filtering
        gender = mycalls[mycalls['library'] == library].loc[:,'gender'].tolist()[0]
        age = mycalls[mycalls['library'] == library].loc[:,'age'].tolist()[0]
        try: #this will only capture values for vafs, depths, and genes if the library is present in df_hopcalls_beataml (the filtered version of df_hopcalls)
            vafs_list = calls[calls['library'] == library].loc[:,'vaf'].tolist()
            depths_list = calls[calls['library'] == library].loc[:,'depth'].tolist()
            genes_list = calls[calls['library'] == library].loc[:,'gene'].tolist()
            ages_list = calls[calls['library'] == library].loc[:,'age'].tolist()
        except: # and creates a blank list if not
            vafs_list = []
            depths_list = []
            genes_list = []
            ages_list = []
        ages.append(age)
        ages_pts_having_calls.append(ages_list)
        vafs.append(vafs_list)
        depths.append(depths_list)
        genes.append(genes_list)
        libs.append(library)
        call_counts.append(len(vafs_list))
        lib_data[library] = [len(vafs_list), age, gender]

        mdn_vafs_ages_zip = zip(mds['ages'], mds['vafs'])
        ages_set = set(mds['ages'])
        mdn_vafs_ages = {}
        for age in ages_set:
            mdn_vafs_ages[age] = []
        for age_i, vaf_i in mdn_vafs_ages_zip:
            mdn_vafs_ages[age_i].append(vaf_i)
        mean_n = np.mean([len(v) for k, v in mdn_vafs_ages.items()])
        mdn_n = np.median([len(v) for k, v in mdn_vafs_ages.items()])

        cnts_ages_zip = zip(ages, call_counts) ## This uses values from df_hopcalls
        ages_set = set(ages)
        cnts_ages = {}
        for age in ages_set:
            cnts_ages[age] = []
        for age_i, cnt_i in cnts_ages_zip:
            cnts_ages[age_i].append(cnt_i)

    x_data = sorted(ages)
    temp_df = pd.DataFrame({'obj':x_data})
    x_data_unique = unique(temp_df['obj'].to_list())
    age_data_labs = [str(e) for e in x_data_unique[:-1]]
    age_data_labs.append('90+')

    # custom labels will change if there are some ages that have no calls 
    ages_w_calls = [item for sublist in ages_pts_having_calls for item in sublist]
    if (len(unique(ages_w_calls)) != len(age_data_labs)): 
        labswithvafs = [str(e) for e in unique(ages_w_calls)[:-1]]
        labswithvafs.append('90+')

    else: 
        labswithvafs = age_data_labs
    


    return {
    'libs': libs,
    'call_counts': call_counts,
    'vafs': vafs,
    'depths': depths,
    'ages': ages,
    'ages_pts_having_calls': ages_pts_having_calls,
    'gends': gends,
    'gends_pts_having_calls': gends_pts_having_calls,
    'genes': genes,
    'lib_data': lib_data, 
    'mean_n': mean_n, 
    'mdn_n': mdn_n, 
    'mdn_vafs_ages' : mdn_vafs_ages, 
    'cnts_ages' : cnts_ages, 
    'age_data_labs': age_data_labs, 
    'labs_w_vafs' : labswithvafs}

def get_frac_pos(all_data):
    age_list = []
    frac_chip_pos = []
    total_list = []
    for i in set(all_data['ages']):
        tot = 0
        pos = 0
        for lib, data in all_data['lib_data'].items():
            if data[1] == i:
                tot += 1
                if data[0] > 0:
                    pos += 1
        age_list.append(i)
        total_list.append(tot)
        frac_chip_pos.append(pos/tot)
    
    return {'ages': age_list, 'fracpos': frac_chip_pos, 'total': total_list}

# get data grouped by gene 
def get_by_gene(calls): 
    genes_mdns = calls.groupby('gene').median(numeric_only=True).sort_values('gene').reset_index()
    genes_means = calls.groupby('gene').mean(numeric_only=True).sort_values('gene').reset_index()
    mdn_depths = genes_mdns.depth.to_list()
    mdn_vafs = genes_mdns.vaf.to_list()
    mean_depth = genes_means.depth.to_list()
    mean_vaf = genes_means.vaf.to_list()

    genes = []
    call_counts = []
    vafs = []
    ages = []
    gends = []
    for gene in sorted(set(calls['gene'])):
        gender = calls[calls['gene'] == gene].loc[:,'gender'].tolist()[0]
        age = calls[calls['gene'] == gene].loc[:,'age'].tolist()[0]
        vafs_list = calls[calls['gene'] == gene].loc[:,'vaf'].tolist()
        ages.append(age)
        vafs.append(vafs_list)
        genes.append(gene)
        call_counts.append(len(vafs_list))

    return {'mdn_depths' : mdn_depths, 'mean_depth' : mean_depth, 'mdn_vafs' : mdn_vafs, 
            'mean_vaf': mean_vaf, 'ages': ages, 'vafs': vafs, 'genes': genes, 'cts':call_counts}



# main function to perform analysis 

def main(): 
    parser = argparse.ArgumentParser(epilog = textwrap.dedent('''\
                                    additional information: 
                                        If you are starting with a file produced from Snakemake, 
                                        --infile, --metadata, --calls, --dbsnp and --plt
                                        are required. Any that you don't fill in will be filled in by default with config file. 
                                        If you have an annotated file already, 
                                        --annot, --dict and --plt are required. 
                                                              
                                    '''))
    
    parser.add_argument("--infile", "-i", help = "Input File produced from Snakemake")
    parser.add_argument("--metadata", "-m", help = "Patient metadata information")
    parser.add_argument("--dbsnp", "-d", help = "Known sites in DBSNP")
    parser.add_argument("--annot", "-a", help = "Annotated file")
    parser.add_argument("--plt", "-p", help = "Directory where plots go")
    parser.add_argument("--ver", "-v", help = "OPTIONAL Can specify version if multiple runs")
    parser.add_argument("--minvaf", help = "OPTIONAL Minimum variant allele frequency")
    parser.add_argument("--maxvaf", help = "OPTIONAL Maximum variant allele frequency")
    parser.add_argument("--mincount", help = "OPTIONAL Minimum read count")
    parser.add_argument("--config", "-c", help="OPTIONAL Path to JSON configuration file", default = '/home/groups/CEDAR/chaoe/HOP_test/config/config_annot.json')
    parser.add_argument("--calls",
                        metavar="KEY=VALUE",
                        nargs='+',
                        help="Set a number of key-value pairs where name = path to calls "), 
    parser.add_argument("--dict", help = "Path to dictionary of calls")

    o = parser.parse_args() 
    
    if o.config: 
        config = load_config(o.config)
        for key, value in config.items(): 
            if not getattr(o, key): 
                setattr(o, key, value)

        print("Pulling default parameters from", o.config)
         
        if o.annot: 
            setattr(o, 'infile', None)   

        
    
    need_annot = validate_args(o)


    while need_annot == True : 
        print("==== Annotations starting ====")
        cols = ['library','chr','pos','ref','alt','depth','count','vaf','avg_BQ','avg_pos_as_fraction','type','class','gene','func', 'dummy']

        # annotate found mutation calls
        df_hopcalls = pd.read_csv(o.infile, sep='\t', index_col=None, names=cols) # without header
        del df_hopcalls['dummy']
        df_hopcalls['gene_ref_aa_pos'] = df_hopcalls.apply(lambda row: [row.gene+'_'+k.replace('p.','')[:-1] for k in re.split('[:,]', row.func) if 'p.' in k], axis = 1)

        # use metadata to add age and gender 
        metadata = pd.read_csv(o.metadata, sep = '\t')
        id_sex_age_dict = {}
        for id in metadata.surgical_pathology_id:
            gender = metadata[(metadata['surgical_pathology_id'] == id)].gender.to_list()[0]
            age = metadata[(metadata['surgical_pathology_id'] == id)].age.to_list()[0]
            id_sex_age_dict[id]= (gender,age)
        
        ages = []
        gends = []
        for index, row in df_hopcalls.iterrows():
            gender = id_sex_age_dict[row.library.rpartition('-')[0]][0]
            age = id_sex_age_dict[row.library.rpartition('-')[0]][1]
            gends.append(gender) 
            ages.append(age) 
        
        df_hopcalls['age'] = ages
        df_hopcalls['gender'] = gends

        print("Added metadata..")

        # annotate also if appears on known mutant lists
        known_calls = dict(map(lambda s: s.split('='), o.calls))

        mt_dict = {k: reformat(v) for k, v in known_calls.items()} 

        matches = []
        for index, row in df_hopcalls.iterrows():
            mutations = set(row['gene_ref_aa_pos'])
            matched_dicts = []
            for k, v in mt_dict.items(): 
                result = not mutations.isdisjoint(v)
                if result: 
                    matched_dicts.append(k)
            if matched_dicts: 
                matches.append(''.join([x[0] for x in matched_dicts]))
            else : 
                matches.append('none')

        df_hopcalls['matches'] = matches

        calls_dict = {k[0] : k for k, v in mt_dict.items()}

        print("Added if calls are known..")
        print("Known calls provided by ", mt_dict.keys())

        # annotate w/ dbsnp id if exists
        use_cols=[0,1,3]
        col_names = ['CHROM','POS', 'ID']
        df_dbsnp = pd.read_csv(o.dbsnp, sep='\t', usecols=use_cols, index_col=None, names=col_names, low_memory=False)
        
        dbsnp_hits = []
        dbsnp_pos = df_dbsnp.POS.to_list()
        for index, row in df_hopcalls.iterrows():
            hit = 'none'
            position = row.iloc[2]
            if position in dbsnp_pos:
                hit = df_dbsnp[df_dbsnp['POS'] == position].loc[:,'ID'].values[0]
            dbsnp_hits.append(hit)

        df_hopcalls['dbsnp'] = dbsnp_hits

        print("Added if dbsnp id exists...")

        # check if diff version number 
        ver = 1 
        if o.ver: 
            ver = o.ver

        # save annotated file and dictionary for calls 
        out_path = os.path.join(os.path.dirname(o.infile), 'HOP_merged_v' + str(ver) + '_' + datetime.today().strftime('%y%m%d') + '.annotated_readcounts_crossref.csv') 
        df_hopcalls.to_csv(out_path, sep='\t', index=False, header = True)

        print("Annotation file has been generated as ", out_path) 

        dict_path = os.path.join(os.path.dirname(o.infile),  'HOP_merged_v' + str(ver) + '_' + datetime.today().strftime('%y%m%d') + '.callsdict.txt') 

        with open(dict_path, "w") as f : 
            json.dump(calls_dict, f)

        print("Dictionary for the calls have been saved as", dict_path )
        need_annot = False 
    
    print("==== Analysis starting ====")

    if o.infile: 
        df_hopcalls = pd.read_csv(out_path, sep = '\t')
    else : 
        df_hopcalls = pd.read_csv(o.annot, sep = '\t')
        dict_path = o.dict

        with open(dict_path, 'r') as f:
            calls_dict = json.load(f)
    
    # pull all data needed for plotting 
    # pull the calls that are present in our annotated df
    diff_calls = set(list(''.join(df_hopcalls[df_hopcalls['matches'] != 'none']['matches']))) 
    calls_analyzed = {}
    # for each unique calls list -> pull the matching data -> values grouped by library, grouped by gene, chip pos calcs etc. 
    for x in diff_calls: 
        calls = df_hopcalls[df_hopcalls['matches'].str.contains(x, na=False)]
        calls, calls_mds, calls_means = get_mdnmean_data(df_hopcalls, calls, o.minvaf, o.maxvaf, o.mincount)
        print(f'{calls_dict[x]} - Number of total calls after filtering: {len(calls)}')
        matching_calls= get_matching_data(df_hopcalls, calls, calls_mds)
        gene_data = [item for sublist in matching_calls['genes'] for item in sublist]
        color_labels = set(gene_data)
        rgb_values = sns.color_palette("tab20",len(color_labels))
        color_map =  dict(zip(color_labels, rgb_values))
        matching_calls['gene_colors'] = [color_map[gene] for gene in gene_data]
        calls_chip_pos = get_frac_pos(matching_calls)
        calls_by_gene = get_by_gene(calls)
        calls_analyzed[calls_dict[x]] = {'calls' : calls, 'mds' : calls_mds, 'means' : calls_means, 'matching' : matching_calls, 'chip_pos' : calls_chip_pos, 'by_gene' : calls_by_gene} 

    
    # plot time!!!
    for k, v in calls_analyzed.items(): 
        path = os.path.join(o.plt, k)
        if not os.path.exists(path): 
            os.mkdir(path)
            os.mkdir(os.path.join(path, 'by_lib')) 
            os.mkdir(os.path.join(path, 'by_age')) 
            os.mkdir(os.path.join(path, 'by_gene')) 
        
        obj = v['matching']
        mds = v['mds']
        means = v['means']
        df = v['calls']
        chip_pos = v['chip_pos']
        by_gene = v['by_gene']

        # plot the ages of each library
        y_data = sorted(obj['ages'])
        x_data = sort_A_using_B(obj['libs'], obj['ages'])
        fig = plt.figure(figsize =(35, 10))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=1)
        plt.title('Age of each Library')
        plt.ylabel('Age (years)')
        plt.xlabel('Library')
        plt.ylim(0,100)
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path,'by_lib/_Ages')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot call counts by library
        y_data = obj['call_counts']
        x_data = obj['libs']
        fig = plt.figure(figsize =(35, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=1)
        plt.title('Call Count of each Library')
        plt.ylabel('Call Count')
        plt.xlabel('Library')
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path,'by_lib/_CallCts')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


        # Plot VAF dist
        y_data = obj['vafs']
        x_data = obj['libs']
        fig = plt.figure(figsize =(35, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('VAF distribution of each Library')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Library')
        plt.margins(x=0)
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path,'by_lib/_VafDist')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot median depth by library
        y_data = mds['unfil_depths']
        x_data = obj['libs']
        fig = plt.figure(figsize =(35, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=1)
        plt.title('Median Depth of each Library')
        plt.ylabel('Median Depth')
        plt.xlabel('Library')
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path,'by_lib/_MedianDepth')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot mean depth by library
        y_data = means['unfil_depths']
        x_data = obj['libs']
        fig = plt.figure(figsize =(35, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=1)
        plt.title('Mean Depth of each Library')
        plt.ylabel('Mean Depth')
        plt.xlabel('Library')
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path,'by_lib/_MeanDepth')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


        # Plot VAF dist by age
        y_data = sort_A_using_B(obj['vafs'], obj['ages'])
        x_data = sorted(obj['ages'])
        fig = plt.figure(figsize =(35, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('VAF distribution of each Library')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Donor Age')
        plt.margins(x=0)
        ax.set_xticklabels(sorted(x_data), rotation=45)#, ha = 'right')
        plt.show()
        file_path = os.path.join(path, 'by_age/_VafDist')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()



        # Plot all vafs by age (scatter) LOG

        y_data = [item for sublist in obj['vafs'] for item in sublist]
        x_data = [item for sublist in obj['ages_pts_having_calls'] for item in sublist]
        z_data = obj['gene_colors']
        labels = unique(x_data[:-1].append('90+')) 
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        jitter(x_data, y_data, c=z_data)
        handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in color_map.items()]
        ax.legend(title='color', handles=handles, bbox_to_anchor=(1, 1), loc='upper left')
        ax.set_xticks(unique(x_data))
        ax.set_xticklabels(obj['labs_w_vafs'])
        plt.yscale('log')
        plt.title('Individual Mutation VAFs by Age (Log)')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Age')
        file_path = os.path.join(path, 'by_age/_VafsLog')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        ## Plot call count by age
        y_data = sort_A_using_B(obj['call_counts'], obj['ages'])
        x_data = sorted(obj['ages'])
        temp_df = pd.DataFrame({'obj':x_data, 'y':y_data})
        fig = plt.figure(figsize =(10, 5))
        ax = sns.stripplot(x=temp_df['obj'], y=temp_df['y'], color='grey', size=8, alpha=.5,linewidth=1.0,edgecolor="black")
        sns.boxplot(x='obj', y='y', data=temp_df,  fliersize=0, color=".8", linecolor="#137", linewidth=.75, medianprops={"color": "r"})
        plt.title('Call counts by Age')
        plt.ylabel('Number of Mutations Passing Filters')
        plt.xlabel(f'Participant Age (years)')
        ax.set_xticks(range(len(unique(x_data))))
        ax.set_xticklabels(obj['age_data_labs'])
        file_path = os.path.join(path, 'by_age/_CallCts')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        
        # Plot call cts by age (scatter)
        y_data = sort_A_using_B(obj['call_counts'], obj['ages'])
        x_data = sorted(obj['ages'])
        fig = plt.figure(figsize =(15, 5))
        plt.title('Age vs Mutations Called')
        plt.ylabel('Mutations Called')
        plt.xlabel('Age')
        #sns.swarmplot(x_data, y_data)
        ax = sns.stripplot(x=x_data, y=y_data, jitter=0.3)
        ax.set_xticks(range(len(unique(x_data))))
        ax.set_xticklabels(obj['age_data_labs'])
        sns.despine()
        file_path = os.path.join(path, 'by_age/_CallCtsScatter')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


        ## Plot call count by age
        y_data = obj['cnts_ages'].values()
        x_data = obj['cnts_ages'].keys()
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('Call counts by Age')
        plt.ylabel('Mutations Called')
        plt.xlabel(f'Library (Median data points for each age = {obj['mdn_n']:.1f})')
        plt.margins(x=0)
        ax.set_xticklabels(sorted(x_data))
        file_path = os.path.join(path, 'by_age/_CallCts2')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


    

        # Plot all vafs by age (scatter)
        y_data = [item for sublist in obj['vafs'] for item in sublist]
        x_data = [item for sublist in obj['ages_pts_having_calls'] for item in sublist]
        z_data = obj['gene_colors']
        fig = plt.figure(figsize =(15, 5))
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        jitter(x_data, y_data, c=z_data)
        handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in color_map.items()]
        ax.legend(title='color', handles=handles, bbox_to_anchor=(1, 1), loc='upper left')
        ax.set_xticks(unique(x_data))
        ax.set_xticklabels(obj['labs_w_vafs'])
        plt.title('Individual Mutation VAFs by Age')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Age')
        file_path = os.path.join(path, 'by_age/_Vafs')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot VAF dist by age
        y_data = obj['mdn_vafs_ages'].values()
        x_data = obj['mdn_vafs_ages'].keys()
    
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('VAF distribution by Age')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel(f'Library (Median data points for each age = {obj['mdn_n']:.1f})')
        plt.margins(x=0)
        ax.set_xticklabels(obj['labs_w_vafs'])
        file_path = os.path.join(path, 'by_age/_VafDist2')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        ## Plot fraction of CHIP-postive donors by age
        y_data = chip_pos['fracpos']
        x_data = chip_pos['ages']
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.scatter(x_data, y_data)
        for i,x in enumerate(x_data):
            ax.text(x-.2, 0.9, chip_pos['total'][i], c='r')
        ax.text(83, 0.95, "Number of Donors:", c='r')
        ax.set_ylim(0,1)
        ax.set_xticks(unique(x_data))
        ax.set_xticklabels(obj['age_data_labs'])
        plt.title('Fraction of CHIP-positive donors by Age')
        plt.ylabel('Fraction of CHIP-Positive Donors')
        plt.xlabel('Donor Age')
        file_path = os.path.join(path, 'by_age/_ChipPos')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

    
        # Plot mean depth by gene
        y_data = by_gene['mean_depth']
        x_data = by_gene['genes']
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=0.8)
        plt.title('Mean Depth of each Gene')
        plt.ylabel('Mean Depth')
        plt.xlabel('Gene')
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path, 'by_gene/_MeanDepth')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

       # Plot call counts by Gene
        y_data = by_gene['cts']
        x_data = by_gene['genes']
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.bar(x_data, y_data, width=0.8)
        plt.title('Call Count of each Gene')
        plt.ylabel('Call Count')
        plt.xlabel('Gene')
        plt.margins(x=0)
        ax.set_xticks(range(len(x_data)))
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path, 'by_gene/_CallCts')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot VAF dist
        y_data = by_gene['vafs']
        x_data = by_gene['genes']
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('VAF distribution of each Gene')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Gene')
        plt.margins(x=0)
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path, 'by_gene/_VAFDist')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot VAF dist LOG
        y_data = by_gene['vafs']
        x_data = by_gene['genes']
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.boxplot(y_data)
        plt.setp(bp['fliers'], color='red', marker='.')
        plt.title('VAF distribution of each Gene')
        plt.ylabel('VAF (<0.45)')
        plt.xlabel('Gene')
        plt.margins(x=0)
        plt.yscale('log')
        ax.set_xticklabels(x_data, rotation=45, ha = 'right')
        file_path = os.path.join(path, 'by_gene/_VAFDistLOG')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


# Plot all vafs by depth (scatter) color coded by gene ID
        y_data = [item for sublist in obj['vafs']for item in sublist]
        x_data = [item for sublist in obj['depths'] for item in sublist]
        z_data = obj['gene_colors']
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.scatter(x_data, y_data, c=z_data)
        handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in color_map.items()]
        ax.legend(title='color', handles=handles, bbox_to_anchor=(1, 1)) #, loc='upper left')
        plt.title('Depth vs VAF')
        plt.ylabel('VAF')
        plt.xlabel('Depth')
        plt.yscale('log')
        file_path = os.path.join(path, '_ScatterVafsbyDepth_GeneColored')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

    

        # Plot call counts by mean depth (scatter)
        y_data = obj['call_counts']
        x_data = mds['unfil_depths']
        corr, _ = pearsonr(y_data, x_data)
        print(k + ' - Pearsons correlation: %.3f' % corr)

        fig = plt.figure(figsize =(7, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.scatter(x_data, y_data)
        plt.title('Mean Depth vs Number of mutations called')
        plt.ylabel('Call Count')
        plt.xlabel('Mean Depth')
        coef = np.polyfit(x_data,y_data,1)
        poly1d_fn = np.poly1d(coef) 
        plt.plot(x_data, poly1d_fn(x_data), '--r')
        corr, _ = pearsonr(y_data, x_data)
        ax.text(1400, .35, f"Pearson's Corr: {corr:.3f}", fontsize=14)#, color='r')
        file_path = os.path.join(path,'_CellCtsbyMeanDepth')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot median vafs by median depth (scatter)
        y_data = mds['vafs']
        x_data = mds['depths']
        fig = plt.figure(figsize =(15, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bp = ax.scatter(x_data, y_data)
        plt.title('Median Depth vs Median VAF')
        plt.ylabel('Median VAF')
        plt.xlabel('Median Depth')
        plt.yscale('log')
        file_path = os.path.join(path, '_MeanVafsbyMeanDepth')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


        # Plot Distribution of Supporting Reads for all Mutations
        x_data = [i for i in df.loc[:,'count'].tolist() if i > 5]
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.hist(x_data,
                bins = 100)
        plt.title('Distribution of Supporting Reads for all Mutations (min 6 supp. reads)')
        plt.ylabel('Count')
        plt.xlabel('Number of supporting reads')
        file_path = os.path.join(path, '_ReadDistMin6')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        x_data = [i for i in df.loc[:,'count'].tolist() if i > 9]
        fig = plt.figure(figsize =(8, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.hist(x_data,
                bins = 100)
        plt.title('Distribution of Supporting Reads for all Mutations (min 10 supp. reads)')
        plt.ylabel('Count')
        plt.xlabel('Number of supporting reads')
        file_path = os.path.join(path, '_ReadDistMin10')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        # Plot distribution of ages across libraries 
        age_counts = Counter(obj['ages'])
        x_data = list(age_counts.keys())
        y_data = list(age_counts.values())
        fig = plt.figure(figsize =(10, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bars = plt.bar(x_data, y_data, color='pink')
        plt.title(f'Distribution of Ages Across Libraries ({len(obj['libs'])} Libraries)' )
        plt.xlabel('Ages')
        plt.ylabel('Number of Libraries')
        ax.set_xticks(unique(x_data))
        ax.set_xticklabels(obj['labs_w_vafs'])
        ax.bar_label(bars)
        file_path = os.path.join(path,'_AgeDistAcrossLibs')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()


        # Plot distribution of call cts across libraries 

        calls = Counter(obj['call_counts'])
        x_data = list(calls.keys())
        y_data = list(calls.values())
        fig = plt.figure(figsize =(10, 5))
        ax = fig.add_axes([0, 0, 1, 1])
        bars = plt.bar(x_data, y_data, color='pink')
        plt.title(f'Distribution of Calls Across Libraries ({len(obj['libs'])} Libraries)' )
        plt.xlabel('Calls')
        plt.ylabel('Number of Libraries')

        ax.bar_label(bars)
        file_path = os.path.join(path,'_CtDistAcrossLibs')
        plt.savefig(file_path, bbox_inches = 'tight')
        plt.close()

        print("Plots generated for", k)







        
main() 