import pandas as pd
import Bio.KEGG.REST as kegg
from collections import OrderedDict
import multiprocessing as mp
import os.path
from functools import partial
import time

io_dir = 'C:\\Users\\Samuel\\Desktop\\postnovo_jpr_materials\\soil_interpretation\\functional_annotation'
input_f = 'unique_ko_entries.tsv'
cores = 4

ko_ec_f = 'ko_ec.tsv'
ec_map_f = 'ec_map.tsv'

def main():

    #Map KO to EC
    df = pd.read_csv(os.path.join(io_dir, input_f), sep='\t', header=None)
    all_kos = []
    for i in df.iloc[:, 0].tolist():
        if ',' in i:
            all_kos += i.split(',')
        else:
            all_kos.append(i)
    unique_kos = sorted(list(set(all_kos)))

    chunk_size = 100
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

    ##Singlethreaded
    #ko_info = []
    #for ko_chunk in ko_chunks:
        #ko_info += kegg_list_worker(ko_chunk)

    #Multithreaded
    pool = mp.Pool(cores)
    ko_info = pool.map(kegg_list_worker, ko_chunks)
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
    pool = mp.Pool(cores)
    map_info = pool.map(kegg_link_worker, ['ec:' + ec for ec in unique_ecs])
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

    return

def kegg_list_worker(query):

    while True:
        try:
            entry = kegg.kegg_list(query).readlines()
            break
        except:
            print('Retrying: ' + query)
            #Biopython says it allows 3 queries per second
            time.sleep(cores * 1/3)

    return entry

def kegg_link_worker(query, db='path'):

    print(query)
    while True:
        try:
            entry = kegg.kegg_link(db, query).readlines()
            break
        except:
            print('Retrying: ' + query)
            #Biopython says it allows 3 queries per second
            time.sleep(cores * 1/3)

    return entry

if __name__ == '__main__':
    main()