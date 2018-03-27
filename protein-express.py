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

script_dir = os.path.dirname(__file__)

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
    '''
    Main method
    '''

    args = get_args()
    global prodigal_bin, blastp_bin, makeblastdb_bin, combined_output_dir
    prodigal_bin = os.path.realpath(args.prodigal)
    blastp_bin = os.path.realpath(args.blastp)
    makeblastdb_bin = os.path.realpath(args.makeblastdb)
    if args.combined_output_dir == None:
        combined_output_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 'search_output'
        )
        try:
            os.mkdir(combined_output_dir)
        except FileExistsError:
            pass
    else:
        combined_output_dir = os.path.realpath(args.combined_output_dir)

    blast_dbs = []
    for pop_bin_fasta in args.population_bins:
        blast_dbs.append(make_blast_db(pop_bin_fasta))

    for proteomic_dataset_dir in args.proteomic_dataset_dirs:
        search_blast_dbs(blast_dbs, proteomic_dataset_dir)

def search_blast_dbs(blast_dbs, proteomic_dataset_dir):
    '''
    Search ORFs of reported peptides against population bin databases 
    '''

    proteome_name = os.path.normpath(os.path.basename(proteomic_dataset_dir))
    query_path = os.path.join(proteomic_dataset_dir, proteome_name + '.blastp_queries.faa')

    if not os.path.exists(query_path):
        reported_peptide_table = os.path.join(proteomic_dataset_dir, 'reported_df.tsv')
        reported_peptide_df = pd.read_csv(reported_peptide_table, sep='\t', header=0)
        scans = reported_peptide_df['scan_list'].tolist()
        # Peptides may have been identified in more than one reference file
        best_predicts_from = reported_peptide_df['best_predicts_from'].tolist()
        best_predicts_from = [s.split(',') for s in best_predicts_from]
        also_contains_predicts_from = reported_peptide_df['also_contains_predicts_from'].tolist()
        also_contains_predicts_from = [
            [''] if pd.isnull(s) else s.split(',') for s in also_contains_predicts_from
        ]
        predicts_from = []
        for i in range(len(best_predicts_from)):
            l = []
            predicts_from.append(l)

            m = best_predicts_from[i]
            if m == ['']:
                pass
            elif '' in m:
                l += m.remove('')
            else:
                l += m

            m = also_contains_predicts_from[i]
            if m == ['']:
                pass
            elif '' in m:
                l += m.remove('')
            else:
                l += m

        # The ORFs containing peptides are in fasta files in the dir for the proteomic dataset
        fasta_files = glob(
            os.path.join(proteomic_dataset_dir, proteome_name + '.*.reads.fasta')
            )
        fasta_files += glob(
            os.path.join(proteomic_dataset_dir, proteome_name + '.*.DBGraphPep2Pro.fasta')
            )
        # Record separately the fasta seq identifiers and seqs from the files
        fasta_identifiers = OrderedDict()
        fasta_seqs = OrderedDict()
        for fasta_file in fasta_files:
            ref_name = os.path.basename(fasta_file).replace(proteome_name + '.', '').replace('.fasta', '')
            with open(fasta_file) as handle:
                lines = handle.readlines()
            fasta_identifiers[ref_name] = [line.rstrip().lstrip('>').split(' ')[0] for line in lines[::2]]
            fasta_seqs[ref_name] = [line.rstrip() for line in lines[1::2]]
            assert len(fasta_identifiers[ref_name]) == len(fasta_seqs[ref_name]), \
            fasta_file + ' does not have an even number of lines'

        # Record the scans and protein matches from each file of db search results
        db_search_results_scans = OrderedDict().fromkeys(fasta_identifiers)
        db_search_results_sources = OrderedDict().fromkeys(fasta_identifiers)
        for ref_name in fasta_identifiers:
            if '.reads' in ref_name:
                db_search_results = proteome_name + '.' + ref_name + \
                '.graph2pro_derep.fgs.tryptic.derep.0.01.tsv'
            elif '.DBGraphPep2Pro' in ref_name:
                db_search_results = proteome_name + '.' + ref_name + '.fixedKR.0.01.tsv'
            db_search_results = os.path.join(proteomic_dataset_dir, db_search_results)
            db_search_results_df = pd.read_csv(db_search_results, sep='\t', header=0)
            db_search_results_scans[ref_name] = db_search_results_df['ScanNum'].tolist()
            l = []
            db_search_results_sources[ref_name] = l
            for protein_names in db_search_results_df['Protein'].tolist():
                l.append([
                    protein_name[:protein_name.index('(pre=')]
                    for protein_name in protein_names.split(';')
                ])

        # Record the ORFs/ORF sources matching each scan in the reported seqs
        first_scans = [int(scans_list_str.split(',')[0]) for scans_list_str in scans]
        orfs = OrderedDict([(first_scan, []) for first_scan in first_scans])
        orf_sources = OrderedDict([(first_scan, []) for first_scan in first_scans])
        for i, scan in enumerate(first_scans):
            # List of non-redundant ORFs matching the scan
            l = orfs[scan]
            m = orf_sources[scan]
            best_predict_sources = best_predicts_from[i]
            for ref_name in best_predict_sources:
                try:
                    j = db_search_results_scans[ref_name].index(scan)
                    for protein_name in db_search_results_sources[ref_name][j]:
                        try:
                            k = fasta_identifiers[ref_name].index(protein_name)
                            orf = fasta_seqs[ref_name][k]
                            try:
                                x = l.index(orf)
                                m[x].append(ref_name)
                            except ValueError:
                                l.append(orf)
                                m.append([ref_name])
                        except ValueError:
                            pass
                except (KeyError, ValueError):
                    pass

        blastp_fasta = []
        for first_scan in first_scans:
            scan_orfs = orfs[first_scan]
            if scan_orfs:
                for j, scan_orf in enumerate(scan_orfs):
                    identifier = '>' + str(first_scan) + '.' + str(j) + '\n'
                    blastp_fasta.append(identifier)
                    blastp_fasta.append(scan_orf + '\n')

        with open(query_path, 'w') as handle:
            for line in blastp_fasta:
                handle.write(line)

    else:
        print(query_path, 'already exists', flush=True)
    
    for blast_db in blast_dbs:
        search_dir = os.path.join(proteomic_dataset_dir, proteome_name + '_pop_bin_search')
        try:
            os.mkdir(search_dir)
        except FileExistsError:
            pass
        out_path = os.path.join(
            search_dir, proteome_name + '.' + os.path.basename(blast_db) + '.blastp_hits.out'
        )
        if not os.path.exists(out_path):
            subprocess.call([
                blastp_bin, 
                '-db', blast_db, 
                '-query', query_path, 
                '-out', out_path, 
                '-evalue', '0.01', 
                '-outfmt', '6', 
                '-comp_based_stats', '0'
            ])
        else:
            print(out_path, 'already exists', flush=True)

        # Load and filter BLAST table
        blast_table = pd.read_csv(out_path, sep='\t', names=blast_table_hdrs, dtype={'qseqid': str})
        blast_table = blast_table.groupby('qseqid', as_index=False).first()
        blast_table = blast_table[blast_table['qend'] - blast_table['length'] == 0]
        blast_table = blast_table[blast_table['evalue'] <= 0.01]
        blast_table = blast_table[blast_table['pident'] >= 95]

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

    return

def make_blast_db(pop_bin_fasta):
    '''
    Make blast database from proteins in the population bin
    '''

    pop_bin_dir = os.path.dirname(pop_bin_fasta)
    pop_bin_name = os.path.splitext(os.path.basename(pop_bin_fasta))[0]
    gene_coords = os.path.join(pop_bin_dir, pop_bin_name + '.gene_coords.gbk')
    proteins = os.path.join(pop_bin_dir, pop_bin_name + '.faa')
    # Run Prodigal to predict genes from the population bin fasta
    if not os.path.exists(proteins):
        subprocess.call([
            prodigal_bin, 
            '-i', pop_bin_fasta, 
            '-o', gene_coords, 
            '-a', proteins
        ])
    else:
        print('prodigal output', proteins, 'already exists', flush=True)

    # Make a blast database from the protein predictions
    blast_db_dir = os.path.join(pop_bin_dir, pop_bin_name + '_blast_db')
    blast_db = os.path.join(blast_db_dir, pop_bin_name)
    try:
        os.mkdir(blast_db_dir)
        subprocess.call([
            makeblastdb_bin, 
            '-dbtype', 'prot', 
            '-in', proteins, 
            '-out', blast_db, 
            '-hash_index'
        ])         
    except FileExistsError:
        print('blast database', blast_db, 'already exists', flush=True)
    
    return blast_db

def get_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--proteomic_dataset_dirs', 
        nargs='+', 
        help=(
            'These directories must be named the same as the dataset and contain: ' + 
            '1. Reported peptides table output by Postnovo' + 
            '2. Filtered database search results ' + 
            '3. Fasta files of the ORFs matched to peptides'
        )
    )
    parser.add_argument(
        '--population_bins', 
        nargs='+', 
        help='Population bin fasta files'
    )
    parser.add_argument(
        '--combined_output_dir', 
        help='Output directory for combined searches against a population bin'
    )
    parser.add_argument(
        '--prodigal', 
        help='prodigal bin'
    )
    parser.add_argument(
        '--makeblastdb', 
        help='makeblastdb bin'
    )
    parser.add_argument(
        '--blastp', 
        help='blastp bin'
    )
    # parser.add_argument(
    #     '--search_result', 
    #     help='Output of this script into which further output will be written'
    # )

    return parser.parse_args()

if __name__ == '__main__':
    main()