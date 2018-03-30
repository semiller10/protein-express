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
from functools import partial
from glob import glob
import multiprocessing as mp
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
bin_table_hdrs = [
    'qfile', 
    'scan', 
    'seqnum', 
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
    'bitscore', 
    'protein', 
    'cog', 
    'descrip'
]

def main():

    args = get_args()
    global out_dir
    out_dir = args.out

    run_prodigal(args.bin_dir, args.threads)

    blast_db_fps = []
    for bin_basename in os.listdir(args.bin_dir):
        blast_db_fps.append(make_blast_db(os.path.join(args.bin_dir, bin_basename)))

    query_fasta_fps = []
    for prot_dir in args.prot_dirs:
        prot_name = os.path.normpath(os.path.basename(prot_dir))
        query_fasta_dir = os.path.join(out_dir, prot_name + '.bin_search')
        query_fasta_fp = os.path.join(query_fasta_dir, prot_name + '.blastp_queries.faa')
        if os.path.exists(query_fasta_fp):
            query_fasta_fps.append(query_fasta_fp)
            print('Query path', query_fasta_fp, 'already exists', flush=True)
        else:
            query_fasta_fps.append(make_query_fasta(prot_dir))

    bin_table_fps = []
    for blast_db_fp in blast_db_fps:
        bin_name = os.path.basename(blast_db_fp)
        bin_table_fp = os.path.join(out_dir, bin_name + '.blast_out.txt')
        bin_table_fps.append(bin_table_fp)
        for i, query_fasta_fp in enumerate(query_fasta_fps):
            prot_name = os.path.basename(query_fasta_fp).replace('.blastp_queries.faa', '')
            search_dir = os.path.join(out_dir, prot_name + '.bin_search')
            blast_table_fp = os.path.join(
                search_dir, prot_name + '.' + os.path.basename(blast_db_fp) + '.blastp_hits.out'
            )
            if os.path.exists(blast_table_fp):
                print(blast_table_fp, 'already exists', flush=True)
            else:
                run_blastp(blast_db_fp, query_fasta_fp, blast_table_fp, args.threads)
            postnovo_table_fp = os.path.join(args.prot_dirs[i], 'reported_df.tsv')
            parse_blast_table(prot_name, blast_table_fp, blast_db_fp, postnovo_table_fp)

    systematize_annot(bin_table_fps)
    compare_bins(bin_table_fps)

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
        '-s', 
        '--state', 
        help='Table relating each proteome name to a state'
    )
    parser.add_argument(
        '-o', 
        '--out', 
        help='Output directory'
    )
    parser.add_argument(
        '-t', 
        '--threads', 
        type=int, 
        default=1, 
        help='Number of threads'
    )

    args = parser.parse_args()

    return args

def run_prodigal(bin_dir, threads):

    bin_basenames = os.listdir(bin_dir)
    bin_fps = [os.path.join(bin_dir, bin_basename) for bin_basename in bin_basenames]

    partial_prodigal_worker = partial(
        prodigal_worker, out_dir=out_dir
    )
    mp_pool = mp.Pool(threads)
    mp_pool.map(partial_prodigal_worker, bin_fps)
    mp_pool.close()
    mp_pool.join()

    return

def prodigal_worker(bin_fp, out_dir):

    bin_dir = os.path.dirname(bin_fp)
    bin_name = os.path.splitext(os.path.basename(bin_fp))[0]
    gene_coords_fp = os.path.join(out_dir, bin_name + '.gene_coords.gbk')
    proteins_fp = os.path.join(out_dir, bin_name + '.faa')

    if not os.path.exists(proteins_fp):
        fnull = open(os.devnull, 'w')
        subprocess.call([
            'prodigal', 
            '-i', bin_fp, 
            '-o', gene_coords_fp, 
            '-a', proteins_fp
        ], stdout=fnull, stderr=fnull)
        fnull.close()
    else:
        print('prodigal output', proteins_fp, 'already exists', flush=True)        

    return

def make_blast_db(bin_fp):
    '''
    Make blast database from proteins in fasta
    '''

    bin_dir = os.path.dirname(bin_fp)
    bin_name = os.path.splitext(os.path.basename(bin_fp))[0]

    # Make a blast database from the proteins
    blast_db_dir = os.path.join(out_dir, bin_name + '.blast_db')
    blast_db_fp = os.path.join(blast_db_dir, bin_name)
    proteins_fp = os.path.join(out_dir, bin_name + '.faa')
    try:
        os.mkdir(blast_db_dir)
        subprocess.call([
            'makeblastdb', 
            '-dbtype', 'prot', 
            '-in', proteins_fp, 
            '-out', blast_db_fp, 
            '-hash_index'
        ])         
    except FileExistsError:
        print('blast database', blast_db_fp, 'already exists', flush=True)

    return blast_db_fp    

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
        [''] if pd.isnull(s) else s.split(',') for s in redun_predict_origins
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
                seq_id = '>' + str(first_scan) + '.' + str(i) + '\n'
                query_fasta.append(seq_id)
                query_fasta.append(orf + '\n')

    search_dir = os.path.join(out_dir, prot_name + '.bin_search')
    subprocess.call(['mkdir', '-p', search_dir])
    query_fasta_fp = os.path.join(search_dir, prot_name + '.blastp_queries.faa')
    with open(query_fasta_fp, 'w') as handle:
        for line in query_fasta:
            handle.write(line)
            
    return query_fasta_fp

def run_blastp(blast_db_fp, query_fasta_fp, blast_table_fp, num_threads):

    print(
        'Aligning', os.path.basename(query_fasta_fp).replace('.blastp_queries.faa', ''), 
        'to', os.path.basename(blast_db_fp), 
        flush=True
    )

    tmp_dir = os.path.join(os.path.dirname(query_fasta_fp), 'tmp')
    query_name = os.path.splitext(os.path.basename(query_fasta_fp))[0]
    subprocess.call(['mkdir', tmp_dir])

    with open(query_fasta_fp) as handle:
        fasta = handle.readlines()

    split_size = len(fasta) // num_threads + ((len(fasta) // num_threads) % 2)
    file_num = 0
    handle = open(os.path.join(tmp_dir, query_name + '.1.faa'), 'w')
    split_fasta_fps = []
    for i, line in enumerate(fasta):
        if i % split_size == 0 and file_num < num_threads:
            file_num += 1
            split_fasta_fp = os.path.join(tmp_dir, query_name + '.' + str(file_num) + '.faa')
            split_fasta_fps.append(split_fasta_fp)
            handle.close()
            handle = open(split_fasta_fp, 'w')
        handle.write(line)
    handle.close()

    partial_blastp_worker = partial(
        blastp_worker, blast_db_fp=blast_db_fp
    )
    mp_pool = mp.Pool(num_threads)
    mp_pool.map(partial_blastp_worker, split_fasta_fps)
    mp_pool.close()
    mp_pool.join()

    blast_df = pd.DataFrame()
    for split_fasta_fp in split_fasta_fps:
        split_blast_table_fp = os.path.splitext(split_fasta_fp)[0] + '.out'
        try:
            split_blast_df = pd.read_csv(split_blast_table_fp, sep='\t', header=None)
            blast_df = pd.concat([blast_df, split_blast_df])
        except pd.io.common.EmptyDataError:
            pass
    blast_df.to_csv(blast_table_fp, sep='\t', index=False, header=False)

    subprocess.call(['rm', '-r', tmp_dir])

    return

def blastp_worker(query_fasta_fp, blast_db_fp):

    split_blast_table_fp = os.path.splitext(query_fasta_fp)[0] + '.out'

    subprocess.call([
        'blastp', 
        '-db', blast_db_fp, 
        '-query', query_fasta_fp, 
        '-out', split_blast_table_fp, 
        '-evalue', '0.01', 
        '-outfmt', '6', 
        '-comp_based_stats', '0'
    ])

    return

def parse_blast_table(prot_name, out_fp, blast_db_fp, postnovo_table_fp):
    '''
    Add BLAST table to merged table for all searches against bin
    '''

    blast_df = pd.read_csv(out_fp, sep='\t', names=blast_table_hdrs, dtype={'qseqid': str})
    blast_df = blast_df[blast_df['length'] >= 9]
    blast_df = blast_df[blast_df['evalue'] <= 0.01]    
    blast_df['scan'] = blast_df['qseqid'].apply(lambda s: int(s.split('.')[0]))
    blast_df['seqnum'] = blast_df['qseqid'].apply(lambda s: s.split('.')[1])
    blast_df.drop('qseqid', axis=1, inplace=True)
    blast_df = blast_df[blast_df.groupby('scan')['evalue'].transform(min) == blast_df['evalue']]
    blast_df = blast_df.groupby('scan', as_index=False).first()
    blast_df.sort_values('scan', inplace=True)

    postnovo_df = pd.read_csv(postnovo_table_fp, sep='\t', header=0)
    postnovo_df['scan'] = postnovo_df['scan_list'].apply(lambda s: int(s.split(',')[0]))
    postnovo_df.sort_values('scan', inplace=True)
    postnovo_df.set_index('scan', inplace=True)
    postnovo_df = postnovo_df.loc[blast_df['scan'].tolist()].reset_index()
    postnovo_df = postnovo_df[['scan', 'predicted name', 'cog cat', 'eggnog hmm desc']]
    postnovo_df.rename(
        columns={'predicted name': 'protein', 'cog cat': 'cog', 'eggnog hmm desc': 'descrip'}, 
        inplace=True
    )
    assert len(blast_df) == len(postnovo_df)
    blast_df = blast_df.merge(postnovo_df, on='scan')
    blast_df['scan'] = blast_df['scan'].apply(str)

    bin_name = os.path.basename(blast_db_fp)
    bin_table_fp = os.path.join(out_dir, bin_name + '.blast_out.txt')
    try:
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)
    except FileNotFoundError:
        bin_df = pd.DataFrame(columns=bin_table_hdrs)

    # Remove any prior entries from the proteomic dataset under consideration
    if prot_name in bin_df['qfile'].tolist():
        bin_df = bin_df[bin_df['qfile'] != prot_name]
    blast_df['qfile'] = prot_name
    bin_df = pd.concat([bin_df, blast_df], ignore_index=True)
    bin_df = bin_df[bin_table_hdrs]
    bin_df.to_csv(bin_table_fp, sep='\t', index=False)

    return bin_table_fp

def systematize_annot(bin_table_fps):
    '''
    Ensure that every protein name maps to the same eggNOG description and COG
    '''

    protein_descrip_counts = dict()
    protein_cog_counts = dict()
    descrip_cog_counts = dict()
    for bin_table_fp in bin_table_fps:
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)[['protein', 'descrip', 'cog']]
        bin_df.fillna('')
        proteins = bin_df['protein'].tolist()
        descrips = bin_df['descrip'].tolist()
        cogs = bin_df['cog'].tolist()
        del(bin_df)

        for i, protein in enumerate(proteins):
            descrip = descrips[i]
            cog = cogs[i]

            if protein != '':
                try:
                    d = protein_descrip_counts[protein]
                    try:
                        d[descrip] += 1
                    except KeyError:
                        d[descrip] = 1
                except KeyError:
                    d = protein_descrip_counts[protein] = dict()
                    d[descrip] = 1

                try:
                    d = protein_cog_counts[protein]
                    try:
                        d[cog] += 1
                    except KeyError:
                        d[cog] = 1
                except KeyError:
                    d = protein_cog_counts[protein] = dict()
                    d[cog] = 1

            try:
                d = descrip_cog_counts[descrip]
                try:
                    d[cog] += 1
                except KeyError:
                    d[cog] = 1
            except KeyError:
                d = descrip_cog_counts[descrip] = dict()
                d[cog] = 1

    protein_descrip = dict()
    protein_cog = dict()
    descrip_cog = dict()
    for protein, descrip_counts in protein_descrip_counts.items():
        max_count = 0
        best_descrip = ''
        for descrip, count in descrip_counts.items():
            if count > max_count:
                best_descrip = descrip
        protein_descrip[protein] = best_descrip
    del(protein_descrip_counts)
    for protein, cog_counts in protein_cog_counts.items():
        max_count = 0
        best_cog = ''
        for cog, count in cog_counts.items():
            if count > max_count:
                best_cog = cog
        protein_cog[protein] = best_cog
    del(protein_cog_counts)
    for descrip, cog_counts in descrip_cog_counts.items():
        max_count = 0
        best_cog = ''
        for cog, count in cog_counts.items():
            if count > max_count:
                best_cog = cog
        descrip_cog[descrip] = best_cog
    del(descrip_cog_counts)

    for bin_table_fp in bin_table_fps:
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)
        bin_df['protein'] = bin_df['protein'].fillna('')
        bin_df['descrip'] = bin_df['descrip'].fillna('')
        bin_df['cog'] = bin_df['cog'].fillna('')
        proteins = bin_df['protein'].tolist()
        old_descrips = bin_df['descrip'].tolist()
        new_descrips = []
        new_cogs = []
        for i, protein in enumerate(proteins):
            if protein == '':
                descrip = old_descrips[i]
                new_descrips.append(descrip)
                new_cogs.append(descrip_cog[descrip])
            else:
                new_descrips.append(protein_descrip[protein])
                new_cogs.append(protein_cog[protein])
        bin_df['descrip'] = new_descrips
        bin_df['cog'] = new_cogs
        bin_df.to_csv(bin_table_fp, sep='\t', index=False)

    return

def compare_bins(bin_table_fps):

    compar_df = pd.DataFrame(columns=['protein', 'descrip', 'cog'])
    for bin_table_fp in bin_table_fps:
        print(bin_table_fp, flush=True)
        bin_name = os.path.basename(bin_table_fp).replace('.blast_out.txt', '')
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)[
            ['protein', 'descrip', 'cog', 'bitscore']
        ]
        protein_gb = bin_df[pd.notnull(bin_df['protein'])].groupby('protein')
        protein_df = protein_gb['bitscore'].agg(['mean', 'min', 'max', 'count'])
        protein_df['protein'] = [protein for protein, _ in protein_gb]
        protein_df['descrip'] = protein_gb['descrip'].agg(lambda d: d.value_counts().index[0])
        protein_df['cog'] = protein_gb['cog'].agg(lambda c: c.value_counts().index[0])
        protein_df.reset_index(inplace=True, drop=True)
        del(protein_gb)
        descrip_gb = bin_df[pd.isnull(bin_df['protein'])].groupby('descrip')
        descrip_df = descrip_gb['bitscore'].agg(['mean', 'min', 'max', 'count'])
        descrip_df['descrip'] = [descrip for descrip, _ in descrip_gb]
        descrip_df['cog'] = descrip_gb['cog'].agg(lambda c: c.value_counts().index[0])
        descrip_df.reset_index(inplace=True, drop=True)
        del(descrip_gb)
        del(bin_df)
        bin_summary_df = pd.concat([protein_df, descrip_df])
        del(protein_df)
        del(descrip_df)
        bin_summary_df = bin_summary_df[
            ['protein', 'descrip', 'cog', 'mean', 'min', 'max', 'count']
        ]
        bin_summary_df['mean'] = bin_summary_df['mean'].round(1)
        bin_summary_df.rename(
            columns={
                'mean': bin_name + '_mean', 
                'min': bin_name + '_min', 
                'max': bin_name + '_max', 
                'count': bin_name + '_count'
            }, 
            inplace=True
        )
        print(bin_table_fp, flush=True)
        compar_df = compar_df.merge(bin_summary_df, how='outer', on=['protein', 'descrip', 'cog'])
        del(bin_summary_df)
    compar_df_fp = os.path.join(out_dir, 'bin_compar.tsv')
    compar_df.sort_values(['protein', 'descrip'], inplace=True)
    compar_df.to_csv(compar_df_fp, sep='\t', index=False)

    return

if __name__ == '__main__':
    main()