#!/usr/bin/env python3

import argparse
import os
import psutil
import sqlite3
import sys

from add_fp_to_db import compose_index_table_name


def sql_for_similarity(fp, mol_field, table, limit=None, radius_morgan=2):
    index_table_name = compose_index_table_name(main_table_name=table, fp=fp, radius_morgan=radius_morgan)
    sql = f"""SELECT 
                    main.smi, 
                    main.id, 
                    bfp_tanimoto(mol_{fp}_bfp(main.{mol_field}, {f'{radius_morgan},' if fp in ['morgan', 'feat_morgan'] else ''} 2048), 
                                 mol_{fp}_bfp(mol_from_smiles(?1), {f'{radius_morgan},' if fp in ['morgan', 'feat_morgan'] else ''} 2048)) as t 
                  FROM 
                    {table} AS main, {index_table_name} AS idx
                  WHERE 
                    main.rowid = idx.id AND
                    idx.id MATCH rdtree_tanimoto(mol_{fp}_bfp(mol_from_smiles(?1), {f'{radius_morgan},' if fp in ['morgan', 'feat_morgan'] else ''} 2048), ?2) 
                  ORDER BY t DESC 
                  {'LIMIT ' + str(limit) if limit is not None else ''}"""
    return sql


def read_smi(fname):
    with open(fname) as f:
        for line in f:
            items = line.strip().split()
            mol_id = items[1] if len(items) > 1 else items[0]
            yield items[0], mol_id


def main():
    parser = argparse.ArgumentParser(description='Similarity search using the selected fingerprints, '
                                                 'which should be previously added to DB.')
    parser.add_argument('-d', '--input_db', metavar='FILENAME', required=True,
                        help='input SQLite DB.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output text file. If omitted output will be printed to STDOUT.')
    parser.add_argument('-q', '--query', metavar='STRING or FNAME', required=True,
                        help='reference SMILES or SMILES file')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where Mol objects are stored. Default: mols.')
    parser.add_argument('-m', '--mol_field', metavar='STRING', default='mol',
                        help='field name where mol objects are stored. Default: mol.')
    parser.add_argument('-f', '--fp', metavar='STRING', default='morgan',
                        choices=['morgan', 'feat_morgan', 'pattern', 'atom_pairs', 'rdkit', 'topological_torsion'],
                        help='fingerprint type to compute. Default: morgan.')
    parser.add_argument('-r', '--radius_morgan', metavar='INTEGER', default=2, type=int,
                        help='radius of Morgan fingerprint. Default: 2.')
    parser.add_argument('-p', '--threshold', metavar='NUMERIC', default=None, type=float,
                        help='Tanimoto similarity threshold. If omitted then top similar compounds will be returned. '
                             'Default: None.')
    parser.add_argument('-l', '--limit', metavar='INTEGER', default=None, type=int,
                        help='maximum number of matches to retrieve. Default: None.')

    args = parser.parse_args()

    con = sqlite3.connect(args.input_db)
    # copy DB to memory
    if psutil.virtual_memory().free > os.path.getsize(args.input_db) * 2:
        dest = sqlite3.connect(':memory:')
        con.backup(dest)
        con.close()
        con = dest

    try:
        con.enable_load_extension(True)
        con.load_extension('chemicalite')
        con.enable_load_extension(False)

        sql = sql_for_similarity(args.fp, args.mol_field, args.table, args.limit, args.radius_morgan)

        if args.output is not None:
            sys.stdout = open(args.output, 'wt')

        if args.threshold is None:
            threshold = 0.7
        else:
            threshold = args.threshold

        sys.stdout.write('\t'.join(['query_smi', 'query_id', 'found_smi', 'found_id', 'similarity']) + '\n')

        if os.path.isfile(args.query):
            for smi, mol_id in read_smi(args.query):
                stop = False
                while not stop:
                    res = con.execute(sql, (smi, threshold)).fetchall()
                    res = [(smi, mol_id) + i for i in res]
                    if res:
                        sys.stdout.write('\n'.join(['\t'.join(map(str, i)) for i in res]) + '\n')
                        sys.stdout.flush()
                        stop = True
                    else:
                        if args.threshold is None:
                            threshold -= 0.1
                        else:
                            stop = True
        else:
            stop = False
            while not stop:
                res = con.execute(sql, (args.query, args.threshold)).fetchall()
                res = [(args.query, args.query) + i for i in res]
                if res:
                    sys.stdout.write('\n'.join(['\t'.join(map(str, i)) for i in res]) + '\n')
                    sys.stdout.flush()
                    stop = True
                else:
                    if args.threshold is None:
                        threshold -= 0.1
                    else:
                        stop = True

    finally:
        con.close()


if __name__ == '__main__':
    main()
