#!/usr/bin/env python3

# Run prodigal on population bin fasta file
# Make blast database from the amino acid fasta output

# Loop through each proteomic dataset
# Loop through each reported peptide (row)
# Loop through each scan in the list of scans that matched the peptide
# Loop through each metagenomic/transcriptomic dataset (best_predict_from and also_contains_predicts_from)
# Search the dataset's filtered database search results for the scan
# Loop through each fasta title
# Recover the ORF from the filtered fasta file
# Add headers/ORFs to a new fasta file
# Headers are in the format: >ScanX.SeqNumY, e.g., >32454.2
# Perform blastp search of the new fasta file against population bin database

# Loop through the blast results
# If the length of the alignment is >9 amino acids, 
# and the e value is <0.01,
# retrieve the scan (query ID) and percent identity
# If at least one peptide matches the population bin, proceed
# Make a file to record results for reported peptides, 
# unless the file exists
# Add a column for the population bin
# Merge the blast results into the table by scan

import argparse
from collections import OrderedDict
from glob import glob
import os
import os.path
import pandas as pd
import subprocess
import sys

blast_table_hdrs = [
    'qseqid', 
    'sseqid', 
    'pident', 
    'length', 
    'mismatch', 
    'gapopen', 
    'qstart', 
    'qend', 
    'sstart', 
    'send', 
    'evalue', 
    'bitscore'
]

def main():

    args = get_args()
    global out_dir = args.out

    blast_db_fps = []
    for fasta_fp in os.listdir(args.bin_dir):
        blast_db_fps.append(make_blast_db(fasta_fp))

    query_fasta_fps = []
    for prot_dir in args.prot_dir:
        prot_name = os.path.normpath(os.path.basename(prot_dir))
        query_fasta_fp = os.path.join(args.out, prot_name + '.blastp_queries.faa')
        if os.path.exists(query_fasta_fp):
            query_fasta_fps.append(query_fasta_fp)
            print('Query path', query_fasta_fp, 'already exists', flush=True)
        else:
            query_fasta_fps.append(make_query_fasta(prot_dir))

    for blast_db_fp in blast_db_fps:
        for query_fasta_fp in query_fasta_fps:
            prot_name = os.path.basename(query_fasta_fp).remove('.blastp_queries.faa')
            search_dir = os.path.join(out_dir, prot_name + '.bin_search')
            out_fp = os.path.join(
                search_dir, prot_name + '.' + os.path.basename(blast_db_fp) + '.blastp_hits.out'
            )
            if os.path.exists(out_fp):
                print(out_fp, 'already exists', flush=True)
            else:
                subprocess.call([
                    blastp_bin, 
                    '-db', blast_db_fp, 
                    '-query', query_fasta_fp, 
                    '-out', out_fp, 
                    '-evalue', '0.01', 
                    '-outfmt', '6', 
                    '-comp_based_stats', '0'
                ])
            parse_blast_table(out_fp, blast_db_fp)

            # Add BLAST table to file of merged tables for the pop bin
            # Remove any prior entries from the proteome under consideration
            pop_bin_name = os.path.basename(blast_db)
            merged_table_path = os.path.join(combined_output_dir, pop_bin_name + '.blast_output.txt')
            merged_table_hdrs = ['qfile'] + blast_table_hdrs[:-1] + blast_table_hdrs[-1:]
            try:
                merged_table = pd.read_csv(merged_table_path, sep='\t', header=0, dtype={'qseqid': str})
            except FileNotFoundError:
                merged_table = pd.DataFrame(columns=merged_table_hdrs)
            if proteome_name in merged_table['qfile'].tolist():
                merged_table = merged_table[merged_table['qfile'] != proteome_name]
            blast_table['qfile'] = proteome_name
            merged_table = pd.concat([merged_table, blast_table], ignore_index=True)
            merged_table = merged_table[merged_table_hdrs]
            merged_table.to_csv(merged_table_path, sep='\t', index=False)    

def parse_blast_table(out_fp, blast_db_fp):

    blast_df = pd.read_csv(out_fp, sep='\t', names=blast_table_hdrs, dtype={'qseqid': str})
    blast_df = blast_df.groupby('qseqid', as_index=False).first()
    blast_df = blast_df[blast_df['length'] >= 9]
    blast_df = blast_df[blast_df['evalue'] <= 0.01]

    bin_name = os.path.basename(blast_db_fp)
    bin_table_fp = os.path.join(out_dir, bin_name + '.blast_out.txt')
    try:
        bin_table_df = pd.read_csv(bin_table_fp, sep='\t', header=0, dtype={'qseqid': str})
    except FileNotFoundError:
        bin_table_df = pd.DataFrame(columns=['qfile'] + blast_table_hdrs)

    return

def make_query_fasta(prot_dir):
    '''
    Recover PSM protein sequences for search against bins
    '''

    postnovo_table_fp = os.path.join(prot_dir, 'reported_df.tsv')
    reported_peptide_df = pd.read_csv(postnovo_table_fp, sep='\t', header=0)
    scans = reported_peptide_df['scan_list'].tolist()
    best_predict_origins = reported_peptide_df['best_predicts_from'].tolist()
    best_predict_origins = [s.split(',') for s in best_predict_origins]
    # Peptides may match a subsequence of a longer ORF from another dataset
    redun_predict_origins = reported_peptide_df['also_contains_predicts_from'].tolist()
    redun_predict_origins = [
        [''] if pd.isnull(s) else s.split(',') for s in redun_predict_origins]
    ]

    predict_origins = []
    # Remove empty strings from origin lists
    for i, best_predict_origin in enumerate(best_predict_origins):
        l = []
        predict_origins.append(l)

        m = best_predict_origin
        if m == ['']:
            pass
        elif '' in m:
            l += m.remove('')
        else:
            l += m

        m = redun_predict_origins[i]
        if m == ['']:
            pass
        elif '' in m:
            l += m.remove('')
        else:
            l += m

    prot_name = os.path.normpath(os.path.basename(prot_dir))
    fasta_fps = glob(
        os.path.join(prot_dir, prot_name + '.*.reads.fasta')
    )
    fasta_fps += glob(
        os.path.join(prot_dir, prot_name + '.*.DBGraphPep2Pro.fasta')
    )
    # Record the seqs and ids from the origin fastas
    fasta_ids = OrderedDict()
    fasta_seqs = OrderedDict()
    for fasta_fp in fasta_fps:
        ref_name = os.path.basename(fasta_fp).replace(prot_name + '.', '').replace('.fasta', '')
        with open(fasta_fp) as handle:
            lines = handle.readlines()
        fasta_ids[ref_name] = [line.lstrip('>').split(' ')[0] for line in lines[::2]]
        fasta_seqs[ref_name] = [line.rstrip() for line in lines[1::2]]
        assert len(fasta_ids[ref_name]) == len(fasta_seqs[ref_name]), \
        fasta_fp + ' does not have an even number of lines'

    # Record the scans and seq matches (hits) from each file of db search results
    db_search_scans = OrderedDict().fromkeys(fasta_ids)
    db_search_hits = OrderedDict().fromkeys(fasta_ids)
    for ref_name in fasta_ids:
        if '.reads' in ref_name:
            db_search_basename = prot_name + '.' + ref_name + \
            '.graph2pro_derep.fgs.tryptic.derep.0.01.tsv'
        elif '.DBGraphPep2Pro' in ref_name:
            db_search_basename = prot_name + '.' + ref_name + '.fixedKR.0.01.tsv'
        db_search_fp = os.path.join(prot_dir, db_search_basename)
        db_search_df = pd.read_csv(db_search_fp, sep='\t', header=0)
        db_search_scans[ref_name] = db_search_df['ScanNum'].tolist()
        l = []
        db_search_hits[ref_name] = l
        for hits_str in db_search_df['Protein'].tolist():
            l.append([
                hit_str[:hit_str.index('(pre=')] for hit_str in hits_str.split(';')
            ])

    first_scans = [int(scans_str.split(',')[0]) for scans_str in scans]
    orfs = OrderedDict([(first_scan, []) for first_scan in first_scans])
    orf_origins = OrderedDict([(first_scan, []) for first_scan in first_scans])
    for i, scan in enumerate(first_scans):
        # List of non-redundant ORFs matching the scan
        orfs_scan = orfs[scan]
        orf_origins_scan = orf_origins[scan]
        best_predict_origins_scan = best_predict_origins[i]
        for ref_name in best_predict_origins_scan:
            # The spectrum (scan) may not have matched the ref db
            try:
                j = db_search_scans[ref_name].index(scan)
                for hit in db_search_hits[ref_name][j]:
                    try:
                        k = fasta_ids[ref_name].index(hit)
                        orf = fasta_seqs[ref_name][k]
                        # The ORF may have already been recorded from a different ref db
                        try:
                            x = orfs_scan.index(orf)
                            orf_origins_scan[x].append(ref_name)
                        except ValueError:
                            orfs_scan.append(orf)
                            orf_origins_scan.append([ref_name])
                    except ValueError:
                        pass
            except (KeyError, ValueError):
                pass

    query_fasta = []
    for first_scan in first_scans:
        orfs_scan = orfs[first_scan]
        if orfs_scan:
            for i, orf in enumerate(orfs_scan):
                seq_id = '>' + str(first_scan) + '.' + str(j) + '\n'
                query_fasta.append(seq_id)
                query_fasta.append(orf + '\n')

    search_dir = os.path.join(out_dir, prot_name + '.bin_search')
    subprocess.call(['mkdir', '-p', search_dir])
    query_fasta_fp = os.path.join(out_dir, search_dir + '.blastp_queries.faa')
    with open(query_fasta_fp, 'w') as handle:
        for line in query_fasta:
            handle.write(line)
            
    return query_fasta_fp

def make_blast_db(fasta_fp):
    '''
    Make blast database from proteins in fasta
    '''

    fasta_dir = os.path.dirname(fasta_fp)
    fasta_name = os.path.splitext(os.path.basename(fasta_fp))[0]
    gene_coords_fp = os.path.join(fasta_dir, fasta_name + '.gene_coords.gbk')
    proteins_fp = os.path.join(fasta_dir, fasta_name + '.faa')
    # Run Prodigal to predict genes
    if not os.path.exists(proteins):
        subprocess.call([
            'prodigal', 
            '-i', fasta_fp, 
            '-o', gene_coords_fp, 
            '-a', proteins_fp
        ])
    else:
        print('prodigal output', proteins_fp, 'already exists', flush=True)

    # Make a blast database from the proteins
    blast_db_dir = os.path.join(out_dir, fasta_name + '.blast_db')
    blast_db_fp = os.path.join(out_dir, fasta_name)
    try:
        os.mkdir(blast_db_dir)
        subprocess.call([
            'makeblastdb', 
            '-dbtype', 'prot', 
            '-in', proteins_f, 
            '-out', blast_db_fp, 
            '-hash_index'
        ])         
    except FileExistsError:
        print('blast database', blast_db_fp, 'already exists', flush=True)
    
    return blast_db_fp

def get_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', 
        '--prot_dirs', 
        nargs='+', 
        help='List of directories for each proteomic dataset'
    )
    parser.add_argument(
        '-b', 
        '--bin_dir', 
        help='Directory exclusively containing bin fastas'
    )
    parser.add_argument(
        '-o', 
        '--out', 
        help='Output directory'
    )

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()