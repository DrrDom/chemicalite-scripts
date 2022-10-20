#!/usr/bin/env python3

import argparse
import pandas as pd
import sqlite3
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import islice
from operator import itemgetter
from sim_search import sql_for_similarity


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def take(n, iterable):
    return list(islice(iterable, n))


def get_similarity(con, fp, mol_field, table, query, threshold, limit, radius_morgan):
        sql = sql_for_similarity(fp=fp, mol_field=mol_field, table=table, limit=limit, radius_morgan=radius_morgan)
        res = con.execute(sql, (query, threshold)).fetchall()
        return res


def calc_sim_for_smiles(smiles, db_name, fp, mol_field, table, threshold, limit, radius_morgan):
    with sqlite3.connect(db_name) as con, sqlite3.connect(':memory:') as dest:
        con.backup(dest)
        dest.enable_load_extension(True)
        dest.load_extension('chemicalite')
        dest.enable_load_extension(False)
        while threshold >= 0:
            all_res = []
            for mol_id, smi in smiles:
                res = get_similarity(dest, fp, mol_field, table, smi, threshold=threshold, limit=limit, radius_morgan=radius_morgan)
                res = [(smi, mol_id) + i for i in res]
                all_res.extend(res)
            df = pd.DataFrame(all_res, columns=['query_smi', 'query_id', 'found_smi', 'found_id', 'similarity'])
            df = df.sort_values('found_id', ascending=False).groupby('found_smi').head(1)
            if len(df) < limit:
                threshold -= 0.05
            else:
                return df.head(limit)
        return None


def main():
    parser = argparse.ArgumentParser(description='Similarity search using the selected fingerprints, '
                                                 'which should be previously added to DB. This script uses all '
                                                 'supplied reference compounds and returns the overall number of '
                                                 'compounds specified in the limit argument across '
                                                 'all given queries together. The threshold will be gradually reduced '
                                                 'starting from the specified value to find the given number of '
                                                 'compounds.')
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
    parser.add_argument('-f', '--fp', metavar='STRING', default='morgan',
                        choices=['morgan', 'feat_morgan', 'pattern', 'atom_pairs', 'rdkit', 'topological_torsion'],
                        help='fingerprint type to compute. Default: morgan.')
    parser.add_argument('-r', '--radius_morgan', metavar='INTEGER', default=2, type=int,
                        help='radius of Morgan fingerprint. Default: 2.')
    parser.add_argument('-p', '--threshold', metavar='NUMERIC', default=0.7, type=float,
                        help='Tanimoto similarity threshold. Default: 0.7.')
    parser.add_argument('-l', '--limit', metavar='INTEGER', required=True, type=int,
                        help='overall maximum number of matches to retrieve across all queries together.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus. This argument has no effect (not implemented).')

    args = parser.parse_args()

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

    df = calc_sim_for_smiles(smiles=list(zip(mol_ids, smiles)),
                             db_name=args.input_db,
                             fp=args.fp,
                             mol_field=args.mol_field,
                             table=args.table,
                             threshold=args.threshold,
                             limit=args.limit,
                             radius_morgan=args.radius_morgan)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == '__main__':
    main()

