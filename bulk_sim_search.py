#!/usr/bin/env python3

import argparse
import sqlite3
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import islice
from sim_search import sql_for_similarity


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def take(n, iterable):
    return list(islice(iterable, n))


def get_similarity(con, fp, mol_field, table, query, limit):
    threshold = 0.7
    while threshold >= 0:
        sql = sql_for_similarity(fp=fp, mol_field=mol_field, table=table, limit=limit)
        res = con.execute(sql, (query, threshold)).fetchall()
        if res:
            return res
        threshold -= 0.1


def calc_sim_for_smiles(smiles, db_name, fp, mol_field, table, limit):
    with sqlite3.connect(db_name) as con, sqlite3.connect(':memory:') as dest:
        con.backup(dest)
        dest.enable_load_extension(True)
        dest.load_extension('chemicalite')
        dest.enable_load_extension(False)
        all_res = []
        for mol_id, smi in smiles:
            res = get_similarity(dest, fp, mol_field, table, smi, limit=limit)
            res = [(smi, mol_id) + i for i in res]
            all_res.extend(res)
    return all_res


def main():
    parser = argparse.ArgumentParser(description='Bulk similarity search using the selected fingerprints, '
                                                 'which should be previously added to DB.')
    parser.add_argument('-d', '--input_db', metavar='FILENAME', required=True,
                        help='input SQLite DB.')
    parser.add_argument('-i', '--input_smiles', metavar='FILENAME', required=True,
                        help='input smiles.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output text file.')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where Mol objects are stored. Default: mols.')
    parser.add_argument('-m', '--mol_field', metavar='STRING', default='mol',
                        help='field name where mol objects are stored. Default: mol.')
    parser.add_argument('-f', '--fp', metavar='STRING', default='morgan', choices=['morgan', 'pattern', 'atom_pairs'],
                        help='fingerprint type to compute. Default: morgan.')
    parser.add_argument('-p', '--threshold', metavar='NUMERIC', default=0.7, type=float,
                        help='Tanimoto similarity threshold. Default: 0.7.')
    parser.add_argument('-l', '--limit', metavar='INTEGER', default=None, type=int,
                        help='maximum number of matches to retrieve. Default: None.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')

    args = parser.parse_args()

    p = Pool(args.ncpu)

    smiles = []
    mol_ids = []
    with open(args.input_smiles) as f:
        for line in f:
            items = line.strip().split()
            smiles.append(items[0])
            if len(items) > 1:
                mol_ids.append(items[1])
            else:
                mol_ids.append(items[0])
    chunked = iter(partial(take, args.ncpu, zip(mol_ids, smiles)), []) # partial calls take function until output is an empty sheet
                                                                       # by take function islice extracts n elements from zip
                                                                       # with saving a condition about zip

    with open(args.output, 'wt') as f:
        f.write('\t'.join(['query_smi', 'query_id', 'found_smi', 'found_id', 'similarity']) + '\n')
        for res in p.map(partial(calc_sim_for_smiles,
                                 db_name=args.input_db,
                                 fp=args.fp,
                                 mol_field=args.mol_field,
                                 table=args.table,
                                 limit=args.limit),
                         chunked):
            for items in res:
                f.write('\t'.join(map(str, items)) + '\n')
                f.flush()


if __name__ == '__main__':
    main()




