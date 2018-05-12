import pandas as pd
import Bio.KEGG.REST as kegg
from collections import OrderedDict
import multiprocessing as mp
import os.path
from functools import partial
import time
import argparse
from glob import glob
import sys
import csv

chunk_size = 100

def main():

    args = parse_args()
    global io_dir, samples_dir, cpus, man_annot_f, ko_ec_f, ec_map_f, map_name_f
    io_dir = args.io_dir
    samples_dir = args.samples_dir
    cpus = args.cpu
    man_annot_f = args.man_annot
    ko_ec_f = args.ko_ec
    ec_map_f = args.ec_map
    map_name_f = args.map_name

    #expand_annotations()
    unique_kos = find_unique_kos()
    ko_ec_dict = ko_to_ec(unique_kos)
    ec_map_dict = ec_to_map(ko_ec_dict)
    map_to_name(ec_map_dict)

    return

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--io_dir', '-d', default='C:\\Users\\Samuel\\Desktop\\postnovo_jpr_materials\\soil_interpretation\\functional_annotation')
    parser.add_argument('--samples_dir', '-s', default='C:\\Users\\Samuel\\Desktop\\postnovo_jpr_materials\\soil_interpretation\\functional_annotation\\sample_reported_dfs')
    parser.add_argument('--cpu', '-c', default=4)
    parser.add_argument('--man_annot', '-a', default='annotations_original_new.tsv')
    parser.add_argument('--ko_ec', '-k', default='ko_ec.tsv')
    parser.add_argument('--ec_map', '-e', default='ec_map.tsv')
    parser.add_argument('--map_name', '-m', default='map_name.tsv')

    return parser.parse_args()

def expand_annotations():

    reannot_dict = OrderedDict()
    with open(os.path.join(io_dir, man_annot_f)) as handle:
        for l in [s.rstrip().split('\t') for s in handle.readlines()]:
            old_name = l[0]
            old_kegg = l[1]
            #Accommodates the following changes to the automatic annotation:
            #1. No name and no KEGG ID were assigned -> New name, no KEGG ID
            #2.                                      -> No name, new KEGG ID
            #3. Name but not KEGG ID were assigned -> Old name, new KEGG ID
            if old_kegg == '':
                if old_name == '':
                    new_name = l[3]
                    new_kegg = l[4]
                    if new_name != '' or new_kegg != '':
                        descrip = l[2]
                        reannot_dict[tuple([old_name, descrip])] = [new_name, new_kegg]
                else:
                    new_kegg = l[4]
                    if new_kegg != '':
                        descrip = l[2]
                        reannot_dict[tuple([old_name, descrip])] = [old_name, new_kegg]

    for sample_dir in glob(os.path.join(samples_dir, '*')):
        reported_df_fp = os.path.join(sample_dir, 'reported_df.tsv')
        reported_df = pd.read_csv(reported_df_fp, sep='\t', header=0)
        old_names = reported_df['predicted name'].fillna('').tolist()
        old_keggs = reported_df['kegg pathways'].fillna('').tolist()
        descrips = reported_df['eggnog hmm desc'].tolist()
        new_names = []
        new_keggs = []
        for old_name, old_kegg, descrip in zip(old_names, old_keggs, descrips):
            if old_kegg == '':
                try:
                    l = reannot_dict[tuple([old_name, descrip])]
                    new_names.append(l[0])
                    new_keggs.append(l[1])
                except KeyError:
                    new_names.append(old_name)
                    new_keggs.append(old_kegg)
            else:
                new_names.append(old_name)
                new_keggs.append(old_kegg)
        reported_df['predicted name'] = new_names
        reported_df['kegg pathways'] = new_keggs

        os.rename(reported_df_fp, os.path.join(sample_dir, 'reported_df.old_annot.tsv'))
        reported_df.to_csv(os.path.join(sample_dir, 'reported_df.tsv'), sep='\t', index=False, quoting=csv.QUOTE_NONE)

    return

def find_unique_kos():

    all_ko_entries = []
    for sample_dir in glob(os.path.join(samples_dir, '*')):
        all_ko_entries = pd.read_csv(
            os.path.join(sample_dir, 'reported_df.tsv'), sep='\t', header=0
        )['kegg pathways'].fillna('').tolist()
    all_kos = []
    for ko_entry in all_ko_entries:
        if ',' in ko_entry:
            all_kos += ko_entry.split(',')
        else:
            all_kos.append(ko_entry)
    unique_kos = sorted(list(set(all_kos)))

    return unique_kos

def ko_to_ec(unique_kos):

    #Map KO to EC
    ko_chunks = []
    ko_chunk = ''
    for i, ko in enumerate(unique_kos):
        if ko_chunk == '':
            ko_chunk = ko
        else:
            ko_chunk += '+' + ko
        if i % chunk_size == chunk_size - 1:
            ko_chunks.append(ko_chunk)
            ko_chunk = ''
        elif i == len(unique_kos) - 1:
            ko_chunks.append(ko_chunk)

    ##Singlethreaded
    #ko_info = []
    #for ko_chunk in ko_chunks:
        #ko_info += kegg_list_worker(ko_chunk)

    #Multithreaded
    pool = mp.Pool(cpus)
    ko_info = pool.map(partial(kegg_list_worker, cpus=cpus), ko_chunks)
    pool.close()
    pool.join()
    #Returns list of lists for each chunk
    ko_info = [i for chunk_info in ko_info for i in chunk_info]

    ko_ec_dict = OrderedDict()
    for info in ko_info:
        #Make sure that the return value has an EC ID at the end of the line in square brackets
        if info[-2:] == ']\n':
            #Queried KO IDs may not have entries in the online database, 
            #so only consider the returned KO IDs
            entry = info.split('\t')
            ko = entry[0].replace('ko:', '')
            ko_ec_dict[ko] = info[info.index('[EC:'):][4: -2].split(' ')

    with open(os.path.join(io_dir, ko_ec_f), 'w') as handle:
        handle.write('KO\tEC\n')
        for ko, ecs in ko_ec_dict.items():
            handle.write(ko + '\t' + ','.join(ecs) + '\n')

    return ko_ec_dict

def kegg_list_worker(query, cpus):

    while True:
        try:
            entry = kegg.kegg_list(query).readlines()
            break
        except:
            print('Retrying: ' + query)
            #Biopython says it allows 3 queries per second
            time.sleep(cpus * 1/3)

    return entry

def ec_to_map(ko_ec_dict):

    #Map EC to Pathway Map
    all_ecs = []
    for ecs in ko_ec_dict.values():
        all_ecs += ecs
    unique_ecs = sorted(list(set(all_ecs)))

    ##Singlethreaded
    #map_info = []
    #for ec in unique_ecs:
    #    map_info += kegg_link_worker('path', 'ec:' + ec).readlines()

    #Multithreaded
    pool = mp.Pool(cpus)
    map_info = pool.map(partial(kegg_link_worker, cpus=cpus), ['ec:' + ec for ec in unique_ecs])
    pool.close()
    pool.join()

    ec_map_dict = OrderedDict()
    for info in map_info:
        ec = info[0].split('\t')[0].replace('ec:', '')
        for entry in info:
            #Returns both Map and EC IDs for each map
            if 'map' in entry:
                map_id = entry.split('\t')[1].replace('path:', '').rstrip()
                try:
                    ec_map_dict[ec].append(map_id)
                except KeyError:
                    ec_map_dict[ec] = [map_id]

    with open(os.path.join(io_dir, ec_map_f), 'w') as handle:
        handle.write('EC\tPathway Map\n')
        for ec, map_ids in ec_map_dict.items():
            handle.write(ec + '\t' + ','.join(map_ids) + '\n')

    return ec_map_dict

def kegg_link_worker(query, cpus, db='path'):

    while True:
        try:
            entry = kegg.kegg_link(db, query).readlines()
            break
        except:
            print('Retrying: ' + query)
            #Biopython says it allows 3 queries per second
            time.sleep(cpus * 1/3)

    return entry

def map_to_name(ec_map_dict):

    #Map Pathway Map ID to Pathway Name
    unique_maps = sorted(list(set(
        [map_id for map_ids in ec_map_dict.values() for map_id in map_ids]
    )))
    map_chunks = []
    map_chunk = ''
    for i, map_id in enumerate(unique_maps):
        if map_chunk == '':
            map_chunk = map_id
        else:
            map_chunk += '+' + map_id
        if i % chunk_size == chunk_size - 1:
            map_chunks.append(map_chunk)
            map_chunk = ''
        elif i == len(unique_maps) - 1:
            map_chunks.append(map_chunk)

    ##Singlethreaded
    #map_info = []
    #for map_chunk in map_chunks:
        #map_info += kegg_list_worker(map_chunk)

    #Multithreaded
    pool = mp.Pool(cpus)
    map_info = pool.map(partial(kegg_list_worker, cpus=cpus), map_chunks)
    pool.close()
    pool.join()
    #Returns list of lists for each chunk
    map_info = [i for chunk_info in map_info for i in chunk_info]

    with open(os.path.join(io_dir, map_name_f), 'w') as handle:
        handle.write('Map ID\tMap Name\n')
        for s in map_info:
            entry = s.split('\t')
            map_id = entry[0].replace('path:', '')
            name = entry[1].rstrip()
            handle.write(map_id + '\t' + name + '\n')

    return

if __name__ == '__main__':
    main()