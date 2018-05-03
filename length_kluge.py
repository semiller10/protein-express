#!/usr/bin/env python3

from collections import OrderedDict
from glob import glob
import numpy as np
import os.path
import pandas as pd
import sys

bin_search_fasta_fp = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/toolik_13_2_1_1/toolik_13_2_1_1.blastp_queries.faa'
sampled_table_fp = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/toolik_13_2_1_1/sampled_table.csv'
emapper_fasta_fps = sorted(glob('/scratch/samuelmiller/12-26-17/postnovo/io/toolik/toolik_13_2_1_1/postnovo_seqs_emapper.*.faa'))
out_dir = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/toolik_13_2_1_1'

def main():

    # Recover first scan IDs of best ORFs used in bin search
    first_scans = ['']
    with open(bin_search_fasta_fp) as handle:
        for i, line in enumerate(handle):
            if i%2 == 0:
                # Header format ex: >36262.0
                first_scan = line.split('.')[0][1:]
                if first_scan != first_scans[-1]:
                    first_scans.append(first_scan)
    first_scans = first_scans[1:]

    # Recover seq numbers corresponding to scan IDs
    sampled_df = pd.read_csv(
        sampled_table_fp, header=0, dtype={'seq_number': str, 'scan_list': str}
    )[['seq_number', 'scan_list']]
    sampled_df['first_scan'] = sampled_df['scan_list'].apply(lambda l: str(sorted(map(int, l.split(',')))[0]))
    sampled_df.set_index('first_scan', inplace=True)
    sampled_df = sampled_df.loc[first_scans, :]
    seq_numbers = sampled_df['seq_number'].tolist()
    seq_number_first_scan_dict = OrderedDict()
    for seq_number, first_scan in zip(seq_numbers, sampled_df.index.tolist()):
        seq_number_first_scan_dict[seq_number] = first_scan
    
    # Load protein subject sequences
    subject_seq_dict = OrderedDict()
    prev_seq_number = ''
    for emapper_fasta_fp in emapper_fasta_fps:
        with open(emapper_fasta_fp) as handle:
            for i, line in enumerate(handle):
                if i%2 == 0:
                    seq_number = line.split('(hit_number)')[0].replace('>(seq_number)', '')
                else:
                    if seq_number == prev_seq_number:
                        subject_seq_dict[seq_number].append(line.rstrip())
                    else:
                        subject_seq_dict[seq_number] = [line.rstrip()]
                        prev_seq_number = seq_number

    avg_protein_lens = []
    std_protein_lens = []
    for seq_number in seq_numbers:
        try:
            protein_lens = list(map(len, subject_seq_dict[seq_number]))
        except KeyError:
            print('Key: ' + str(seq_number))
            sys.exit()
        avg_protein_lens.append(sum(protein_lens))
        std_protein_lens.append(np.std(protein_lens))
    pd.DataFrame(
        OrderedDict([
            ('seq_number', seq_numbers), 
            ('first_scan', [seq_number_first_scan_dict[seq_number] for seq_number in seq_numbers]), 
            ('avg_protein_len', avg_protein_lens), 
            ('std_protein_len', std_protein_lens)
        ])
    ).to_csv(os.path.join(out_dir, 'avg_protein_lens.tsv'), sep='\t', index=False)

    return

def add_len_col():

    reported_df = pd.read_csv(os.path.join(out_dir, 'reported_df.tsv'), sep='\t', header=0)
    reported_df['first_scan'] = reported_df['scan_list'].apply(lambda l: l.split(',')[0])
    lens_df = pd.read_csv(
        os.path.join(out_dir, 'avg_protein_lens.tsv'), 
        dtype={'first_scan': str, 'avg_protein_len': str}, 
        sep='\t', 
        header=0
    )
    lens_df.set_index('first_scan', inplace=True)
    lens_df = lens_df[['avg_protein_len']]
    lens_df.columns = ['protein length']
    lens_df['protein length'] = lens_df['protein length'].apply(int).apply(str)
    lens_df = lens_df.groupby(lens_df.index).first()
    new_header = reported_df.columns.tolist()
    new_header.insert(new_header.index('scan count') + 1, 'protein length')
    new_reported_df = reported_df.merge(lens_df, how='outer', left_on=reported_df['first_scan'], right_index=True)
    new_reported_df = new_reported_df[new_header]
    new_reported_df.to_csv(os.path.join(out_dir, 'reported_df2.tsv'), sep='\t', index=False)

    return

if __name__ == '__main__':
    #main()
    add_len_col()