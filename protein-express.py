#!/usr/bin/env python3

import argparse
import Bio.KEGG.REST as kegg
from collections import OrderedDict
from copy import deepcopy
from functools import partial
from glob import glob
import multiprocessing as mp
import numpy as np
import os
import os.path
import pandas as pd
import pickle as pkl
import seaborn as sb
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
trans_table = dict.fromkeys(
    map(ord, ''.join(['.', '|', '^', '+', '-'] + [str(i) for i in range(10)])), None
    )
ranks = [
    'species', 
    'genus', 
    'family', 
    'order', 
    'class', 
    'phylum', 
    'superkingdom'
    ]
bin_table_hdrs = [
    'qfile', 
    'scan', 
    'seqnum', 
    'peptide', 
    'scans', 
    'speccount', 
    'protlen', 
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
    ] + ranks
compar_table_merge_hdrs = ['protein', 'descrip', 'cog']

sig_cutoff = 27

def main():

    args = get_args()
    global prot_dirs, bin_dir, out_dir, num_threads
    prot_dirs = args.prot_dirs
    bin_dir = args.bin_dir
    out_dir = args.out
    num_threads = args.threads

    global kegg_dir
    if args.kegg_dir != None:
        kegg_dir = args.kegg_dir
    else:
        kegg_dir = None

    prot_state = OrderedDict()
    if args.state != None:
        state_df = pd.read_csv(args.state, sep='\t', header=0)
        for _, row in state_df.iterrows():
            prot = row.iloc[0]
            state = row.iloc[1]
            prot_state[prot] = state
        global all_states
        all_states = list(set(state_df.iloc[:, 1].tolist()))

    # Predict genes from bin contigs
    run_prodigal()

    # Make BLAST+ database from each bin
    blast_db_fps = []
    for bin_basename in os.listdir(bin_dir):
        blast_db_fps.append(make_blast_db(os.path.join(bin_dir, bin_basename)))
    global bin_names
    bin_names = [os.path.basename(blast_db_fp) for blast_db_fp in blast_db_fps]

    # BLAST PSM ORFs against bins
    query_fasta_fps = []
    global prot_names
    prot_names = []
    for prot_dir in args.prot_dirs:
        prot_name = os.path.normpath(os.path.basename(prot_dir))
        prot_names.append(prot_name)
        query_fasta_dir = os.path.join(out_dir, prot_name + '.bin_search')
        query_fasta_fps.append(os.path.join(query_fasta_dir, prot_name + '.blastp_queries.faa'))
        if os.path.exists(query_fasta_fps[-1]):
            print('Query path', query_fasta_fps[-1], 'already exists', flush=True)
        else:
            make_query_fasta(prot_name, prot_dir)

    if prot_state:
        state_pipeline(prot_state, blast_db_fps, query_fasta_fps)
    else:
        stateless_pipeline(blast_db_fps, query_fasta_fps)
        
    if args.group != None:
        score_df = assign_groups(args.group)
        
    return

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
    parser.add_argument(
        '-k', 
        '--kegg_dir', 
        help='Directory containing ko_ec.tsv, ec_map.tsv, map_name.tsv'
    )
    parser.add_argument(
        '-s', 
        '--state', 
        help='Table relating each proteome name to a state'
    )
    parser.add_argument(
        '-g', 
        '--group', 
        help='Table relating protein names or descriptions to groups'
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

def run_prodigal():

    bin_basenames = os.listdir(bin_dir)
    bin_fps = [os.path.join(bin_dir, bin_basename) for bin_basename in bin_basenames]

    mp_pool = mp.Pool(num_threads)
    mp_pool.map(prodigal_worker, bin_fps)
    mp_pool.close()
    mp_pool.join()

    return

def prodigal_worker(bin_fp):

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

def make_query_fasta(prot_name, prot_dir):
    '''
    Recover PSM ORF sequences for search against bins
    '''

    postnovo_table_fp = os.path.join(prot_dir, 'reported_df.tsv')
    reported_peptide_df = pd.read_csv(postnovo_table_fp, sep='\t', header=0)
    scans = reported_peptide_df['scan_list'].tolist()
    # Dataset (e.g., metagenome) origins of unique ORFs containing PSMs
    best_predict_origins = reported_peptide_df['best_predicts_from'].tolist()
    best_predict_origins = [s.split(',') for s in best_predict_origins]
    # Origins of non-unique ORFs with PSMs that are subseqs of unique ORFs
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

    # Get the names of all origin datasets
    fasta_fps = glob(
        os.path.join(prot_dir, prot_name + '.*.reads.fasta')
    )
    fasta_fps += glob(
        os.path.join(prot_dir, prot_name + '.*.DBGraphPep2Pro.fasta')
    )
    # Record the ORF seqs and ids from the origin fastas
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

    # Record the scans, peptides, and ORF seqs from each file of db search results
    db_search_scans = OrderedDict().fromkeys(fasta_ids)
    db_search_peps = OrderedDict().fromkeys(fasta_ids)
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
        db_search_peps[ref_name] = db_search_df['Peptide'].tolist()
        l = []
        db_search_hits[ref_name] = l
        for hits_str in db_search_df['Protein'].tolist():
            l.append([
                hit_str[:hit_str.index('(pre=')] for hit_str in hits_str.split(';')
            ])

    first_scans = [int(scans_str.split(',')[0]) for scans_str in scans]
    # It is possible that a spectrum may have been ID'd as different precursor peptide sequences
    peps = OrderedDict([(first_scan, []) for first_scan in first_scans])
    orfs = OrderedDict([(first_scan, []) for first_scan in first_scans])
    orf_origins = OrderedDict([(first_scan, []) for first_scan in first_scans])
    for i, scan in enumerate(first_scans):
        # List of non-redundant ORFs matching the scan
        peps_scan = peps[scan]
        orfs_scan = orfs[scan]
        orf_origins_scan = orf_origins[scan]
        best_predict_origins_scan = best_predict_origins[i]
        for ref_name in best_predict_origins_scan:
            # The spectrum (scan) may not have matched the ref db
            try:
                j = db_search_scans[ref_name].index(scan)
                db_search_ref_peps = db_search_peps[ref_name]
                for hit in db_search_hits[ref_name][j]:
                    # I don't remember why I included the next exception, and don't think it occurs
                    try:
                        k = fasta_ids[ref_name].index(hit)
                        orf = fasta_seqs[ref_name][k]
                        # The ORF may have already been recorded from a different ref db
                        try:
                            x = orfs_scan.index(orf)
                            orf_origins_scan[x].append(ref_name)
                        except ValueError:
                            # Assume that peptides are the same if ORFs are the same
                            peps_scan.append(db_search_ref_peps[j])
                            orfs_scan.append(orf)
                            orf_origins_scan.append([ref_name])
                    except ValueError:
                        pass
            except (KeyError, ValueError):
                pass

    search_dir = os.path.join(out_dir, prot_name + '.bin_search')
    subprocess.call(['mkdir', '-p', search_dir])

    peps_fp = os.path.join(search_dir, prot_name + '.peps.pkl')
    with open(peps_fp, 'wb') as handle:
        pkl.dump(peps, handle, 2)

    query_fasta = []
    for first_scan in first_scans:
        orfs_scan = orfs[first_scan]
        if orfs_scan:
            for i, orf in enumerate(orfs_scan):
                # Index the ORFs containing the same peptide
                seq_id = '>' + str(first_scan) + '.' + str(i) + '\n'
                query_fasta.append(seq_id)
                query_fasta.append(orf + '\n')

    query_fasta_fp = os.path.join(search_dir, prot_name + '.blastp_queries.faa')
    with open(query_fasta_fp, 'w') as handle:
        for line in query_fasta:
            handle.write(line)
            
    return

def state_pipeline(prot_state, blast_db_fps, query_fasta_fps):

    state_bin_table_fp = OrderedDict([(state, []) for state in all_states])
    bin_table_fps = []
    for i, bin_name in enumerate(bin_names):
        blast_db_fp = blast_db_fps[i]
        bin_table_fps_for_bin = []
        for state in all_states:
            bin_table_fp = os.path.join(out_dir, bin_name + '.' + state + '.blast_out.txt')
            state_bin_table_fp[state].append(bin_table_fp)
            bin_table_fps_for_bin.append(bin_table_fp)
        bin_table_fps += bin_table_fps_for_bin
        for prot_name, query_fasta_fp, prot_dir in zip(prot_names, query_fasta_fps, prot_dirs):
            state = prot_state[prot_name]
            search_dir = os.path.join(out_dir, prot_name + '.bin_search')
            blast_table_fp = os.path.join(
                search_dir, prot_name + '.' + bin_name + '.blastp_hits.out'
            )
            if os.path.exists(blast_table_fp):
                print(blast_table_fp, 'already exists', flush=True)
            else:
                run_blastp(blast_db_fp, query_fasta_fp, blast_table_fp)
            postnovo_table_fp = os.path.join(prot_dir, 'reported_df.tsv')
            peps_fp = os.path.join(search_dir, prot_name + '.peps.pkl')
            ## UNCOMMENT
            #parse_blast_table(
            #    prot_name, 
            #    blast_table_fp, 
            #    blast_db_fp, 
            #    postnovo_table_fp, 
            #    peps_fp, 
            #    os.path.join(out_dir, bin_name + '.' + state + '.blast_out.txt')
            #)
        for bin_table_fp in bin_table_fps_for_bin:
            remove_redun_peps(bin_table_fp)

    ## UNCOMMENT
    #systematize_annot(bin_table_fps)
    state_compar_table_fp = OrderedDict().fromkeys(state_bin_table_fp)
    ## UNCOMMENT
    #for state, bin_table_fps_for_state in state_bin_table_fp.items():
    #    compar_table_basename = 'compar_table.' + state + '.tsv'
    #    state_compar_table_fp[state] = compare_bins(
    #        bin_table_fps_for_state, compar_table_basename
    #    )
    state_compar_table_fp['tussock'] = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/protein-express_out/compar_table.tussock.tsv'
    state_compar_table_fp['intertussock'] = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/protein-express_out/compar_table.intertussock.tsv'
    state_compar_table_fp['shrub'] = '/scratch/samuelmiller/12-26-17/postnovo/io/toolik/protein-express_out/compar_table.shrub.tsv'
    compare_states(state_compar_table_fp)

    #Make tables of NSAF statistics for proteins in KEGG maps
    if kegg_dir != None:
        make_vanted_map_tables(prot_state)

    return

def stateless_pipeline(blast_db_fps, query_fasta_fps):

    bin_table_fps = []
    # Search each bin
    for i, bin_name in enumerate(bin_names):
        blast_db_fp = blast_db_fps[i]
        bin_table_fp = os.path.join(out_dir, bin_name + '.blast_out.txt')
        bin_table_fps.append(bin_table_fp)
        for prot_name, query_fasta_fp, prot_dir in zip(prot_names, query_fasta_fps, prot_dirs):
            search_dir = os.path.join(out_dir, prot_name + '.bin_search')
            blast_table_fp = os.path.join(
                search_dir, prot_name + '.' + bin_name + '.blastp_hits.out'
            )
            if os.path.exists(blast_table_fp):
                print(blast_table_fp, 'already exists', flush=True)
            else:
                run_blastp(blast_db_fp, query_fasta_fp, blast_table_fp)
            postnovo_table_fp = os.path.join(prot_dir, 'reported_df.tsv')
            peps_fp = os.path.join(search_dir, prot_name + '.peps.pkl')
            parse_blast_table(
                prot_name, blast_table_fp, blast_db_fp, postnovo_table_fp, peps_fp, bin_table_fp
            )
        remove_redun_peps(bin_table_fps[-1])

    systematize_annot(bin_table_fps)
    compare_bins(bin_table_fps)

    #Make tables of NSAF statistics for proteins in KEGG maps
    if kegg_dir != None:
        prot_state = OrderedDict()
        for prot_name in prot_names:
            prot_state[prot_name] = prot_name
        global all_states
        all_states = list(set(prot_state.values()))
        make_vanted_map_tables(prot_state)

    return

def vanted_col_setup():

    #Make table formatted for plotting pathway map graphs in Vanted
    vanted_cols = OrderedDict()
    vanted_cols[0] = [''] * 22
    vanted_cols[0][0] = 'VANTED - Input File'
    vanted_cols[0][2] = 'Experiment'
    vanted_cols[0][3] = 'Start of Experiment (Date)'
    vanted_cols[0][4] = 'Remark*'
    vanted_cols[0][5] = 'Experiment Name (ID)'
    vanted_cols[0][6] = 'Coordinator'
    vanted_cols[0][7] = 'Sequence-Name*'
    vanted_cols[0][10] = 'Plants/Genotypes**'
    vanted_cols[0][11] = 'Species'
    vanted_cols[0][12] = 'Variety*'
    vanted_cols[0][13] = 'Genotype'
    vanted_cols[0][14] = 'Growth conditions*'
    vanted_cols[0][15] = 'Treatment*'
    vanted_cols[0][19] = 'Measurements'
    vanted_cols[0][21] = 'Plant/Genotype***'

    vanted_cols[1] = [''] * 22
    vanted_cols[1][3] = '1/1/2018'
    vanted_cols[1][5] = 'Experiment'
    vanted_cols[1][6] = 'Scientist'
    vanted_cols[1][21] = 'Replicate #'

    vanted_cols[2] = [''] * 22
    vanted_cols[2][21] = 'Time*'

    vanted_cols[3] = [''] * 22
    vanted_cols[3][21] = 'Unit (Time)*'

    vanted_cols[4] = [''] * 22
    vanted_cols[4][2] = 'Important Info'
    vanted_cols[4][3] = '- Fields with a * are optional'
    vanted_cols[4][4] = '- Yellow cells allow input'
    vanted_cols[4][5] = '** These cells must contain numbers as 1, 2, 3, ...'
    vanted_cols[4][6] = '*** These cells must correlate to the numbers in **'
    vanted_cols[4][7] = '- The Experiment Name must be unique in the whole database'
    vanted_cols[4][19] = 'Substance'
    vanted_cols[4][20] = 'Meas.-Tool*'
    vanted_cols[4][21] = 'Unit'

    vanted_cols[5] = [''] * 22
    vanted_cols[6] = [''] * 22
    vanted_cols[7] = [''] * 22
    vanted_cols[8] = [''] * 22
    vanted_cols[9] = [''] * 22
    vanted_cols[10] = [''] * 22
    vanted_cols[10][2] = 'Internal Info'
    vanted_cols[10][3] = 'v1.2'

    if len(all_states) > 11:
        for i in range(11, len(all_states)):
            vanted_cols[i] = [''] * 22

    vanted_row_names = OrderedDict()
    vanted_row_names['Plants/Genotypes**'] = 10
    vanted_row_names['Species'] = 11
    vanted_row_names['Genotype'] = 13

    return vanted_cols, vanted_row_names

def make_vanted_map_tables(prot_state_dict):

    vanted_cols, vanted_row_names = vanted_col_setup()

    #Mapping of KO IDs to lists of EC IDs
    ko_ec_dict = OrderedDict()
    with open(os.path.join(kegg_dir, 'ko_ec.tsv')) as handle:
        for entry in [s.rstrip().split('\t') for s in handle.readlines()[1:]]:
            ko_ec_dict[entry[0]] = entry[1].split(',')
    #Mapping of EC IDs to Map IDs
    ec_map_dict = OrderedDict()
    with open(os.path.join(kegg_dir, 'ec_map.tsv')) as handle:
        for entry in [s.rstrip().split('\t') for s in handle.readlines()[1:]]:
            ec_map_dict[entry[0]] = entry[1].split(',')
    #Mapping of Map IDs to Pathway Names
    map_name_dict = OrderedDict()
    with open(os.path.join(kegg_dir, 'map_name.tsv')) as handle:
        for entry in [s.rstrip().split('\t') for s in handle.readlines()[1:]]:
            map_name_dict[entry[0]] = entry[1]

    #Three strategies for assigning SAF to protein functional annotations:
    #Only consider peptides with KEGG KO IDs, grouping results by KEGG KO ID
    #Only consider peptides with gene family names, grouping results by gene family name
    #Consider peptides with KEGG KO IDs, 
    #then consider peptides with gene family names lacking KEGG KO IDs, 
    #grouping by KEGG KO ID or gene family name

    #Mapping of EC ID to score for each proteomic sample
    ec_saf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    sample_ec_nsaf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    state_ec_nsaf_dicts = OrderedDict(
        [(state, OrderedDict()) for state in all_states]
    )
    #Mapping of gene family name to score for each proteomic sample
    gene_saf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    sample_gene_nsaf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    state_gene_nsaf_dicts = OrderedDict(
        [(state, OrderedDict()) for state in all_states]
    )
    #Mapping of gene family name (for peptides without KEGG KO IDs) to score for each proteomic sample
    gene_wout_ec_saf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    sample_gene_wout_ec_nsaf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    state_gene_wout_ec_nsaf_dicts = OrderedDict(
        [(state, OrderedDict()) for state in all_states]
    )
    #Mapping of Map ID to summed EC NSAF values for each state
    sample_map_nsaf_dicts = OrderedDict(
        [(prot_name, OrderedDict()) for prot_name in prot_state_dict]
    )
    state_map_nsaf_dicts = OrderedDict(
        [(state, OrderedDict()) for state in all_states]
    )
    #Mapping of Map ID to EC IDs
    map_ec_dict = OrderedDict()
    #Spectrum count per sample
    sample_count_dict = OrderedDict([(prot_name, 0) for prot_name in prot_state_dict])
    #Spectrum count per state
    state_count_dict = OrderedDict([(state, 0) for state in all_states])
    
    #Add general Vanted input table information for each state
    for i, state in enumerate(all_states):
        vanted_cols[i + 1][vanted_row_names['Plants/Genotypes**']] = i + 1
        vanted_cols[i + 1][vanted_row_names['Species']] = state
        vanted_cols[i + 1][vanted_row_names['Genotype']] = 'Mixed'
        vanted_cols[0].append(i + 1)
        vanted_cols[1].append(1)
        vanted_cols[2].append(0)
        vanted_cols[3].append('day')

    for prot_name, prot_dir in zip(prot_names, prot_dirs):
        # Calculate spectral count statistics for every unique peptide in a sample
        postnovo_table_fp = os.path.join(prot_dir, 'reported_df.tsv')
        postnovo_df = pd.read_csv(postnovo_table_fp, sep='\t', header=0)
        postnovo_df = postnovo_df[pd.notnull(postnovo_df['protein length'])]
        postnovo_df['SAF'] = postnovo_df['scan count'] / postnovo_df['protein length']
        sample_count_dict[prot_name] = postnovo_df['scan count'].sum()
        state_count_dict[prot_state_dict[prot_name]] += sample_count_dict[prot_name]

        #Peptides with KO IDs
        kegg_df = postnovo_df[pd.notnull(postnovo_df['kegg pathways'])]
        #Peptides with gene family name or KO IDs
        gene_kegg_df = postnovo_df[
            pd.notnull(postnovo_df['predicted name']) | 
            pd.notnull(postnovo_df['kegg pathways'])
        ]

        #Record sample SAF values for peptides with gene family name
        gene_saf_dict = gene_saf_dicts[prot_name]
        gene_wout_ec_saf_dict = gene_wout_ec_saf_dicts[prot_name]
        for gene_name, ko_entry, saf in zip(
            gene_kegg_df['predicted name'].fillna('').tolist(), 
            gene_kegg_df['kegg pathways'].fillna('').tolist(), 
            gene_kegg_df['SAF'].fillna('').tolist()
        ):
            if ko_entry == '':
                try:
                    gene_wout_ec_saf_dict[gene_name] += saf
                except KeyError:
                    gene_wout_ec_saf_dict[gene_name] = saf
            try:
                gene_saf_dict[gene_name] += saf
            except KeyError:
                gene_saf_dict[gene_name] = saf

        #Record sample SAF values for peptides with KEGG KO IDs
        ec_saf_dict = ec_saf_dicts[prot_name]
        for ko_entry, saf in zip(kegg_df['kegg pathways'].tolist(), kegg_df['SAF'].tolist()):
            kos = ko_entry.split(',')
            ko_count = len(kos)
            for ko in kos:
                try:
                    for ec in ko_ec_dict[ko]:
                        #When a peptide has multiple KO IDs, 
                        #treat each of these as a separate guess at the true function, 
                        #so scale the SAF by the number of KO IDs for the entry
                        adjusted_saf = saf / ko_count
                        try:
                            ec_saf_dict[ec] += adjusted_saf
                        except KeyError:
                            ec_saf_dict[ec] = adjusted_saf
                #KO ID doesn't map to an EC ID
                except KeyError:
                    pass

    #Two ways of turning SAF into NSAF (determining which peptide count to use as denominator):
    #By sample
    #By class of sample (state) -- this treats samples in a state as a single "metasample"
    for prot_name, state in prot_state_dict.items():

        ec_saf_dict = ec_saf_dicts[prot_name]
        sample_ec_nsaf_dict = sample_ec_nsaf_dicts[prot_name]
        state_ec_nsaf_dict = state_ec_nsaf_dicts[state]
        sample_map_nsaf_dict = sample_map_nsaf_dicts[prot_name]
        state_map_nsaf_dict = state_map_nsaf_dicts[state]
        for ec, saf in ec_saf_dict.items():
            sample_nsaf = saf / sample_count_dict[prot_name]
            state_nsaf = saf / state_count_dict[state]
            sample_ec_nsaf_dict[ec] = sample_nsaf
            try:
                state_ec_nsaf_dict[ec] += state_nsaf
            except KeyError:
                state_ec_nsaf_dict[ec] = state_nsaf

            #Determine the Pathway Maps associated with the EC IDs
            try:
                for map_id in ec_map_dict[ec]:
                    try:
                        map_ec_dict[map_id].append(ec)
                    except KeyError:
                        map_ec_dict[map_id] = [ec]
                    try:
                        sample_map_nsaf_dict[map_id] += sample_nsaf
                    except KeyError:
                        sample_map_nsaf_dict[map_id] = sample_nsaf
                    try:
                        state_map_nsaf_dict[map_id] += state_nsaf
                    except KeyError:
                        state_map_nsaf_dict[map_id] = state_nsaf
            #EC ID doesn't map to a Map ID
            except KeyError:
                pass

        gene_saf_dict = gene_saf_dicts[prot_name]
        sample_gene_nsaf_dict = sample_gene_nsaf_dicts[prot_name]
        state_gene_nsaf_dict = state_gene_nsaf_dicts[state]
        for gene, saf in gene_saf_dict.items():
            sample_gene_nsaf_dict[gene] = saf / sample_count_dict[prot_name]
            try:
                state_gene_nsaf_dict[gene] += saf / state_count_dict[state]
            except KeyError:
                state_gene_nsaf_dict[gene] = saf / state_count_dict[state]

        gene_wout_ec_saf_dict = gene_wout_ec_saf_dicts[prot_name]
        sample_gene_wout_ec_nsaf_dict = sample_gene_wout_ec_nsaf_dicts[prot_name]
        state_gene_wout_ec_nsaf_dict = state_gene_wout_ec_nsaf_dicts[state]
        for gene, saf in gene_wout_ec_saf_dict.items():
            sample_gene_wout_ec_nsaf_dict[gene] = saf / sample_count_dict[prot_name]
            try:
                state_gene_wout_ec_nsaf_dict[gene] += saf / state_count_dict[state]
            except KeyError:
                state_gene_wout_ec_nsaf_dict[gene] = saf / state_count_dict[state]

    #Dereplicate list of sample EC IDs for each Map ID
    for map_id, ecs in map_ec_dict.items():
        map_ec_dict[map_id] = list(set(ecs))

    #Make a Vanted-formatted table for each Map ID, and add info for each EC ID
    map_table_dict = OrderedDict()
    map_nsaf_items = []
    for state_map_nsaf_dict in state_map_nsaf_dicts.values():
        map_nsaf_items += state_map_nsaf_dict.items()
    unique_map_ids = [t[0] for t in sorted(map_nsaf_items, key=lambda t: -t[1])]
    for map_id in unique_map_ids:
        map_table_dict[map_id] = deepcopy(vanted_cols)
        map_vanted_cols = map_table_dict[map_id]
        for i, ec in enumerate(map_ec_dict[map_id]):
            if 5 + i >= len(map_vanted_cols):
                map_vanted_cols[5 + i] = [''] * 22
            map_vanted_cols[5 + i][19] = ec
            map_vanted_cols[5 + i][20] = 'MS'
            map_vanted_cols[5 + i][21] = 'NSAF Score'
            for state in all_states:
                try:
                    nsaf = state_ec_nsaf_dicts[state][ec]
                except KeyError:
                    nsaf = 0
                map_vanted_cols[5 + i].append(str(nsaf))

    #Write the Vanted-formatted tables to Excel spreadsheets
    kegg_map_dir = os.path.join(out_dir, 'vanted_kegg_map_data')
    if not os.path.exists(kegg_map_dir):
        os.mkdir(kegg_map_dir)
    total_map_nsaf_data = []
    for map_id, vanted_map_cols in map_table_dict.items():
        #Equalize column lengths with empty cells
        max_len = max([len(col) for col in vanted_map_cols.values()])
        for i, col in vanted_map_cols.items():
            while len(col) < max_len:
                col.append('')
        vanted_df = pd.DataFrame(vanted_map_cols)
        map_name = map_name_dict[map_id]
        #Replace backslashes in Map Names, e.g., "Glycolysis / Gluconeogenesis"
        vanted_f = map_name.replace('/', '-') + '.xlsx'
        vanted_fp = os.path.join(kegg_map_dir, vanted_f)
        writer = pd.ExcelWriter(vanted_fp, engine='xlsxwriter')
        vanted_df.to_excel(writer, sheet_name='VANTED Input Form', header=False, index=False)
        date_format = writer.book.add_format()
        date_format.set_num_format('dd/mm/yyyy')
        writer.sheets['VANTED Input Form'].write(3, 1, '1/1/2018', date_format)
        writer.save()
        #Write a table of each pathway's total score from all states into the map data directory
        total_map_nsaf = 0
        for state, state_map_nsaf_dict in state_map_nsaf_dicts.items():
            try:
                total_map_nsaf += state_map_nsaf_dict[map_id]
            except KeyError:
                pass
        total_map_nsaf_data.append([map_id, map_name, total_map_nsaf])
    total_map_nsaf_df = pd.DataFrame(
        total_map_nsaf_data, columns=['Pathway Map', 'Pathway Map Name', 'Total Score']
    )
    total_map_nsaf_fp = os.path.join(kegg_map_dir, 'total_map_nsaf.tsv')
    total_map_nsaf_df.to_csv(total_map_nsaf_fp, sep='\t', index=False)

    #Write tables of NSAF values for multivariate analysis
    nsaf_data_dir = os.path.join(out_dir, 'nsaf_data')
    if not os.path.exists(nsaf_data_dir):
        os.mkdir(nsaf_data_dir)
    for prot_name, state in prot_state_dict.items():

        ec_sample_nsaf_f = 'EC_NSAF.samples.tsv'
        ec_sample_nsaf_fp = os.path.join(nsaf_data_dir, ec_sample_nsaf_f)
        additional_df = pd.DataFrame(list(sample_ec_nsaf_dicts[prot_name].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays(
            [[state], [prot_name]], names=['State', 'Sample']
        )
        if os.path.exists(ec_sample_nsaf_fp):
            existing_df = pd.read_csv(ec_sample_nsaf_fp, sep='\t', index_col=0, header=[0, 1])
            merged_df = pd.merge(
                existing_df, additional_df, how='outer', left_index=True, right_index=True
            )   
        else:
            merged_df = additional_df
        merged_df.to_csv(ec_sample_nsaf_fp, sep='\t')

        map_sample_nsaf_f = 'MAP_NSAF.samples.tsv'
        map_sample_nsaf_fp = os.path.join(nsaf_data_dir, map_sample_nsaf_f)
        additional_df = pd.DataFrame(list(sample_map_nsaf_dicts[prot_name].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays(
            [[state], [prot_name]], names=['State', 'Sample']
        )
        if os.path.exists(map_sample_nsaf_fp):
            existing_df = pd.read_csv(map_sample_nsaf_fp, sep='\t', index_col=0, header=[0, 1])
            merged_df = pd.merge(
                existing_df, additional_df, how='outer', left_index=True, right_index=True
            )   
        else:
            merged_df = additional_df
        merged_df.to_csv(map_sample_nsaf_fp, sep='\t')

        gene_sample_nsaf_f = 'GENE_NSAF.samples.tsv'
        gene_sample_nsaf_fp = os.path.join(nsaf_data_dir, gene_sample_nsaf_f)
        additional_df = pd.DataFrame(list(sample_gene_nsaf_dicts[prot_name].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays(
            [[state], [prot_name]], names=['State', 'Sample']
        )
        if os.path.exists(gene_sample_nsaf_fp):
            existing_df = pd.read_csv(gene_sample_nsaf_fp, sep='\t', index_col=0, header=[0, 1])
            merged_df = pd.merge(
                existing_df, additional_df, how='outer', left_index=True, right_index=True
            )
        else:
            merged_df = additional_df
        merged_df.to_csv(gene_sample_nsaf_fp, sep='\t')

        map_plus_gene_sample_nsaf_f = 'MAP_GENE_NSAF.samples.tsv'
        map_plus_gene_sample_nsaf_fp = os.path.join(nsaf_data_dir, map_plus_gene_sample_nsaf_f)
        additional_df = pd.concat([
            pd.DataFrame(list(sample_ec_nsaf_dicts[prot_name].items())).set_index(0), 
            pd.DataFrame(list(sample_gene_wout_ec_nsaf_dicts[prot_name].items())).set_index(0)
        ])
        additional_df.columns = pd.MultiIndex.from_arrays(
            [[state], [prot_name]], names=['State', 'Sample']
        )
        if os.path.exists(map_plus_gene_sample_nsaf_fp):
            existing_df = pd.read_csv(
                map_plus_gene_sample_nsaf_fp, sep='\t', index_col=0, header=[0, 1]
            )
            merged_df = pd.merge(
                existing_df, additional_df, how='outer', left_index=True, right_index=True
            )
        else:
            merged_df = additional_df
        merged_df.to_csv(map_plus_gene_sample_nsaf_fp, sep='\t')

    #Sort columns to ensure that states are grouped
    for fp in [
        ec_sample_nsaf_fp, map_sample_nsaf_fp, gene_sample_nsaf_fp, map_plus_gene_sample_nsaf_fp
    ]:
        df = pd.read_csv(fp, sep='\t', index_col=0, header=[0, 1])
        df.sort_index(axis=1, level=['State', 'Sample'], inplace=True)
        df.to_csv(fp, sep='\t')

    for state in all_states:

        ec_state_nsaf_f = 'EC_NSAF.metasample.tsv'
        ec_state_nsaf_fp = os.path.join(nsaf_data_dir, ec_state_nsaf_f)
        additional_df = pd.DataFrame(list(state_ec_nsaf_dicts[state].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays([[state]], names=['State'])
        if os.path.exists(ec_state_nsaf_fp):
            existing_df = pd.read_csv(ec_state_nsaf_fp, sep='\t', index_col=0, header=0)
            merged_df = pd.merge(existing_df, additional_df, how='outer', left_index=True, right_index=True)   
        else:
            merged_df = additional_df
        merged_df.to_csv(ec_state_nsaf_fp, sep='\t')

        map_state_nsaf_f = 'MAP_NSAF.metasample.tsv'
        map_state_nsaf_fp = os.path.join(nsaf_data_dir, map_state_nsaf_f)
        additional_df = pd.DataFrame(list(state_map_nsaf_dicts[state].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays([[state]], names=['State'])
        if os.path.exists(map_state_nsaf_fp):
            existing_df = pd.read_csv(map_state_nsaf_fp, sep='\t', index_col=0, header=0)
            merged_df = pd.merge(existing_df, additional_df, how='outer', left_index=True, right_index=True)   
        else:
            merged_df = additional_df
        merged_df.to_csv(map_state_nsaf_fp, sep='\t')

        gene_state_nsaf_f = 'GENE_NSAF.metasample.tsv'
        gene_state_nsaf_fp = os.path.join(nsaf_data_dir, gene_state_nsaf_f)
        additional_df = pd.DataFrame(list(state_gene_nsaf_dicts[state].items())).set_index(0)
        additional_df.columns = pd.MultiIndex.from_arrays([[state]], names=['State'])
        if os.path.exists(gene_state_nsaf_fp):
            existing_df = pd.read_csv(gene_state_nsaf_fp, sep='\t', index_col=0, header=0)
            merged_df = pd.merge(existing_df, additional_df, how='outer', left_index=True, right_index=True)   
        else:
            merged_df = additional_df
        merged_df.to_csv(gene_state_nsaf_fp, sep='\t')

        map_plus_gene_state_nsaf_f = 'MAP_GENE_NSAF.metasample.tsv'
        map_plus_gene_state_nsaf_fp = os.path.join(nsaf_data_dir, map_plus_gene_state_nsaf_f)
        additional_df = pd.concat([
            pd.DataFrame(list(state_ec_nsaf_dicts[state].items())).set_index(0), 
            pd.DataFrame(list(state_gene_wout_ec_nsaf_dicts[state].items())).set_index(0)
        ])
        additional_df.columns = pd.MultiIndex.from_arrays([[state]], names=['State'])
        if os.path.exists(map_plus_gene_state_nsaf_fp):
            existing_df = pd.read_csv(map_plus_gene_state_nsaf_fp, sep='\t', index_col=0, header=0)
            merged_df = pd.merge(existing_df, additional_df, how='outer', left_index=True, right_index=True)
        else:
            merged_df = additional_df
        merged_df.to_csv(map_plus_gene_state_nsaf_fp, sep='\t')

    ##Temp
    #with open(os.path.join(kegg_map_dir, 'general_data.txt'), 'w') as handle:
    #    for state in all_states:
    #        handle.write('Count of spectra with KEGG annotation in ' + state + ' samples: ' + str(spec_count_dict[state]) + '\n')
    #    handle.write('\n')
    #    for state in all_states:
    #        handle.write('EC IDs with highest NSAF in ' + state + ' samples\n')
    #        state_ec_score_dict = state_ec_score_dicts[state]
    #        for t in sorted(state_ec_score_dict.items(), key=lambda t: -t[1])[81:100]:
    #            try:
    #                enzyme_name = kegg.kegg_list('ec:' + t[0]).readlines()[0].split('\t')[1].split(';')[0]
    #            except:
    #                enzyme_name = ''
    #            handle.write(t[0] + ' ' + enzyme_name + ' ' + str(t[1]) + '\n')
    #        handle.write('\n')

    return

def run_blastp(blast_db_fp, query_fasta_fp, blast_table_fp):

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

def parse_blast_table(prot_name, out_fp, blast_db_fp, postnovo_table_fp, peps_fp, bin_table_fp):
    '''
    Add BLAST table to merged table for all searches against bin
    '''

    blast_df = pd.read_csv(out_fp, sep='\t', names=blast_table_hdrs, dtype={'qseqid': str})
    #BLAST+ is currently set to report alignments with an evalue <= 0.01
    #blast_df = blast_df[blast_df['evalue'] <= 0.01]
    blast_df['scan'] = blast_df['qseqid'].apply(lambda s: int(s.split('.')[0]))
    blast_df['seqnum'] = blast_df['qseqid'].apply(lambda s: int(s.split('.')[1]))
    blast_df.drop('qseqid', axis=1, inplace=True)
    #Retain the ORF with the lowest evalue for each peptide
    blast_df = blast_df[blast_df.groupby('scan')['evalue'].transform(min) == blast_df['evalue']]
    blast_df = blast_df.groupby('scan', as_index=False).first()
    blast_df.sort_values('scan', inplace=True)

    postnovo_df = pd.read_csv(postnovo_table_fp, dtype={'scan count': str}, sep='\t', header=0)
    # Each peptide is ID'd by the first scan in the scan list
    postnovo_df['scan'] = postnovo_df['scan_list'].apply(lambda s: int(s.split(',')[0]))
    postnovo_df.sort_values('scan', inplace=True)
    postnovo_df.set_index('scan', inplace=True)
    postnovo_df = postnovo_df.loc[blast_df['scan'].tolist()].reset_index()
    postnovo_df = postnovo_df[
        [
            'scan', 
            'scan_list', 
            'scan count', 
            'protein length', 
            'predicted name', 
            'cog cat', 
            'eggnog hmm desc'
        ] + ranks
    ]
    postnovo_df.rename(
        columns={
            'scan_list': 'scans', 
            'scan count': 'speccount', 
            'protein length': 'protlen', 
            'predicted name': 'protein', 
            'cog cat': 'cog', 
            'eggnog hmm desc': 'descrip'
        }, 
        inplace=True
    )
    assert len(blast_df) == len(postnovo_df)
    merged_df = blast_df.merge(postnovo_df, on='scan')
    del(blast_df)

    with open(peps_fp, 'rb') as handle:
        peps = pkl.load(handle)
    orf_peps = []
    seq_nums = merged_df['seqnum'].tolist()
    for scan, seq_num in zip(merged_df['scan'].tolist(), seq_nums):
        orf_peps.append(peps[scan][seq_num].translate(trans_table))
    merged_df['peptide'] = orf_peps

    merged_df['scan'] = merged_df['scan'].apply(str)

    try:
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)
    except FileNotFoundError:
        bin_df = pd.DataFrame(columns=bin_table_hdrs)

    # Remove any prior entries from the proteomic dataset under consideration
    if prot_name in bin_df['qfile'].tolist():
        bin_df = bin_df[bin_df['qfile'] != prot_name]
    merged_df['qfile'] = prot_name
    bin_df = pd.concat([bin_df, merged_df], ignore_index=True)
    bin_df = bin_df[bin_table_hdrs]
    bin_df.to_csv(bin_table_fp, sep='\t', index=False)

    return bin_table_fp

def remove_redun_peps(bin_table_fp):
    
    #When a peptide is found in multiple proteomic datasets, 
    #all instances are compared to the same reference data, 
    #so all instances have the same ORFs, 
    #so all instances have the same bin alignment results.
    bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)
    bin_df = bin_df.groupby('peptide', as_index=False).first()
    bin_df = bin_df[bin_table_hdrs]
    bin_df.sort_values(['qfile', 'scan'], inplace=True)
    bin_df.to_csv(bin_table_fp, sep='\t', index=False)

    return

def systematize_annot(bin_table_fps):
    '''
    Ensure that protein names map to consistent eggNOG descriptions and COGs
    '''

    protein_descrip_counts = dict()
    protein_cog_counts = dict()
    descrip_cog_counts = dict()
    for bin_table_fp in bin_table_fps:
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)[compar_table_merge_hdrs]
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

def compare_bins(bin_table_fps, compar_table_basename='bin_compar.tsv'):

    compar_df = pd.DataFrame(columns=compar_table_merge_hdrs)
    for i, bin_table_fp in enumerate(bin_table_fps):
        bin_name = bin_names[i]
        bin_df = pd.read_csv(bin_table_fp, sep='\t', header=0)[
            ['protein', 'descrip', 'cog', 'bitscore'] + ranks
        ]
        # Group peptide results for bin into protein results
        # First consider peptides with protein name annotations
        protein_gb = bin_df[pd.notnull(bin_df['protein'])].groupby('protein')
        protein_df = protein_gb['bitscore'].agg(['mean', 'min', 'max', 'count'])
        protein_df['protein'] = [protein for protein, _ in protein_gb]
        protein_df['descrip'] = protein_gb['descrip'].agg(lambda d: d.value_counts().index[0])
        protein_df['cog'] = protein_gb['cog'].agg(lambda c: c.value_counts().index[0])
        protein_df.reset_index(inplace=True, drop=True)
        protein_df[ranks] = merge_ranks(protein_gb)
        del(protein_gb)
        # Second consider peptides that only have description string annotations
        bin_df['descrip_lowercase'] = bin_df['descrip'].str.lower()
        descrip_gb = bin_df[pd.isnull(bin_df['protein'])].groupby('descrip_lowercase')
        descrip_df = descrip_gb['bitscore'].agg(['mean', 'min', 'max', 'count'])
        descrip_df['descrip'] = descrip_gb['descrip'].agg(lambda d: d.value_counts().index[0])
        descrip_df['cog'] = descrip_gb['cog'].agg(lambda c: c.value_counts().index[0])
        descrip_df.reset_index(inplace=True, drop=True)
        descrip_df[ranks] = merge_ranks(descrip_gb)
        del(descrip_gb)
        del(bin_df)
        # Summarize the information from each bin in a table, with a block of columns for each bin
        bin_summary_df = pd.concat([protein_df, descrip_df], ignore_index=True)
        del(protein_df)
        del(descrip_df)
        bin_summary_df = bin_summary_df[
            ['protein', 'descrip', 'cog', 'mean', 'min', 'max', 'count'] + ranks
        ]
        bin_summary_df['mean'] = bin_summary_df['mean'].round(1)
        bin_summary_df.rename(
            columns=dict(
                (old_name, bin_name + '_' + old_name) for old_name 
                in ['mean', 'min', 'max', 'count'] + ranks
            ), 
            inplace=True
        )
        compar_df = compar_df.merge(bin_summary_df, how='outer', on=compar_table_merge_hdrs)
        del(bin_summary_df)

    # Find an LCA lineage of the protein from the PSM-ORF data
    ranks_df = pd.DataFrame(columns=compar_table_merge_hdrs + ranks)
    for bin_name in bin_names:
        bin_rank_cols = [bin_name + '_' + rank for rank in ranks]
        bin_ranks_df = compar_df[compar_table_merge_hdrs + bin_rank_cols].copy()
        bin_ranks_df.rename(
            columns=dict(
                (bin_rank_col, ranks[i]) for i, bin_rank_col in enumerate(bin_rank_cols)
            ), 
            inplace=True
        )
        ranks_df = pd.concat([ranks_df, bin_ranks_df], ignore_index=True)
        compar_df.drop(bin_rank_cols, axis=1, inplace=True)
    compar_df['protein'] = compar_df['protein'].fillna('')
    compar_df.sort_values(compar_table_merge_hdrs, inplace=True)
    compar_df.reset_index(drop=True, inplace=True)
    ranks_df['protein'] = ranks_df['protein'].fillna('')
    compar_df[ranks] = merge_ranks(ranks_df.groupby(compar_table_merge_hdrs, as_index=False))

    # Count the number of peptides in each bin that were found for the protein
    compar_df['total_count'] = 0
    for bin_name in bin_names:
        bin_count_hdr = bin_name + '_count'
        compar_df[bin_count_hdr].fillna(0, inplace=True)
        compar_df['total_count'] += compar_df[bin_count_hdr]
        compar_df[bin_count_hdr].replace(0, np.nan, inplace=True)

    compar_df.sort_values(['cog', 'total_count'], ascending=[True, False], inplace=True)
    compar_table_fp = os.path.join(out_dir, compar_table_basename)
    compar_df.to_csv(compar_table_fp, sep='\t', index=False)

    return compar_table_fp

def merge_ranks(gb):

    agg_tax = OrderedDict([(rank, []) for rank in ranks])
    for _, group_df in gb:
        tax_consistency = False
        for rank in ranks:
            if tax_consistency:
                agg_tax[rank].append(group_df[rank].iloc[0])
            else:
                rank_series = group_df[rank]
                if len(rank_series.unique()) == 1:
                    tax = rank_series.iloc[0]
                    if pd.isnull(tax):
                        agg_tax[rank].append('')
                    else:
                        tax_consistency = True
                        agg_tax[rank].append(tax)
                else:
                    agg_tax[rank].append('')
    ranks_df = pd.DataFrame.from_dict(agg_tax)

    return ranks_df

def compare_states(state_compar_table_fp):

    state_compar_table = OrderedDict([
        (state, pd.read_csv(compar_table_fp, sep='\t', header=0)) 
        for state, compar_table_fp in state_compar_table_fp.items()
    ])

    # Merge sample columns from each state's sample comparison table
    num_cols = len(list(state_compar_table.values())[0].columns)
    all_state_df = pd.DataFrame(columns=compar_table_merge_hdrs)
    i = 0
    while True:
        prev_state_suff = '_' + list(state_compar_table.keys())[0]
        col_df = list(state_compar_table.values())[0]
        col = col_df.columns[i + len(compar_table_merge_hdrs)]
        col_df = col_df[compar_table_merge_hdrs + [col]]
        col_df.rename(columns={col: col + prev_state_suff}, inplace=True)
        for state, compar_df in list(state_compar_table.items())[1:]:
            state_suff = '_' + state
            next_col_df = compar_df[compar_table_merge_hdrs + [col]]
            next_col_df.rename(columns={col: col + state_suff}, inplace=True)
            col_df = col_df.merge(
                next_col_df, 
                how='outer', 
                on=compar_table_merge_hdrs, 
                suffixes=(prev_state_suff, state_suff)
            )
            prev_state_suff = state_suff
        all_state_df = all_state_df.merge(
            col_df, how='outer', on=compar_table_merge_hdrs
        )
        i += 1
        if i == num_cols - len(compar_table_merge_hdrs):
            break
            
    for state in all_states:
        all_state_df[state + '_max_mean'] = all_state_df[
            [bin_name + '_mean_' + state for bin_name in bin_names]
        ].max(1)

    total_count_cols = all_state_df[
        ['total_count_' + state for state in all_states]
    ].copy().fillna(0)
    all_state_df['total_count_diff'] = total_count_cols.max(1) - total_count_cols.min(1)

    all_state_df.sort_values(['cog', 'total_count_diff'], ascending=[True, False], inplace=True)
    all_state_table_fp = os.path.join(out_dir, 'state_compar.tsv')
    all_state_df.to_csv(all_state_table_fp, sep='\t', index=False)

    return

def assign_groups(assign_fp):
    '''
    Assign proteins and descriptions to groups
    '''

    assign_df = pd.read_csv(
        assign_fp, sep='\t', header=0, names=['group', 'protein', 'descrip']
        )
    assign_df.fillna('', inplace=True)
    protein_group = OrderedDict([
        (p, g) for p, g 
        in zip(assign_df['protein'].tolist(), assign_df['group'])
        if p != ''
        ])
    descrip_group = OrderedDict([
        (d, g) for d, g 
        in zip(assign_df['descrip'].tolist(), assign_df['group'])
        if d != ''
        ])
    
    for state in all_states:
        fp = os.path.join(out_dir, 'compar_table.' + state + '.tsv')
        compar_df = pd.read_csv(fp, sep='\t', header=0)
        compar_df['protein'] = compar_df['protein'].fillna('')
        compar_df = compar_df[
            (compar_df['protein'].isin(protein_group))
            | ((compar_df['protein'] == '') & (compar_df['descrip'].isin(descrip_group)))
            ]
        protein_df = pd.DataFrame(columns=['group', 'protein', 'descrip'])
        protein_df['protein'] = compar_df['protein']
        protein_df['descrip'] = compar_df['descrip']
        groups = []
        for p, d in zip(compar_df['protein'].tolist(), compar_df['descrip'].tolist()):
            if p == '':
                groups.append(descrip_group[d])
            else:
                groups.append(protein_group[p])
        protein_df['group'] = groups

        for bin in bin_names:
            bin_scores = []
            bin_df = compar_df[[bin + '_mean', bin + '_count']]
            bin_df.fillna(0, inplace=True)
            bin_df[bin + '_mean'] = bin_df[bin + '_mean'] - sig_cutoff
            protein_df[bin] = bin_df[bin + '_mean'] * bin_df[bin + '_count']

        protein_df.sort_values(['group', 'protein'], inplace=True)
        protein_df.to_csv(
            os.path.join(out_dir, 'protein_score.' + state + '.tsv'), sep='\t', index=False
            )

        group_df = protein_df.groupby('group').sum()
        score_df = group_df.div(group_df.max(0), axis=1)
        score_df.to_csv(
            os.path.join(out_dir, 'bin.group_score.' + state + '.tsv'), sep='\t'
            )

    return score_df

if __name__ == '__main__':
    main()