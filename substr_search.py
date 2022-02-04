#!/usr/bin/env python3

import argparse
import sqlite3
import sys


def main():
    parser = argparse.ArgumentParser(description='Substructure search using pattern fingerprints, '
                                                 'which should be previously added to DB.')
    parser.add_argument('-d', '--input_db', metavar='FILENAME', required=True,
                        help='input SQLite DB.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output text file. If omitted output will be printed to STDOUT.')
    parser.add_argument('-q', '--query', metavar='STRING', required=True,
                        help='reference SMILES.')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where Mol objects are stored. Default: mols.')
    parser.add_argument('-m', '--mol_field', metavar='STRING', default='mol',
                        help='field name where mol objects are stored. Default: mol.')
    parser.add_argument('-l', '--limit', metavar='INTEGER', default=None, type=int,
                        help='maximum number of matches to retrieve. Default: None.')

    args = parser.parse_args()
    with sqlite3.connect(args.input_db) as con:

        con.enable_load_extension(True)
        con.load_extension('chemicalite')
        con.enable_load_extension(False)

        sql = f"""SELECT {args.table}.smi, {args.table}.id FROM {args.table}, {args.table}_pattern_idx AS idx WHERE
                  {args.table}.rowid = idx.id AND
                  mol_is_substruct({args.table}.{args.mol_field}, mol_from_smiles(?1)) AND
                  idx.id MATCH rdtree_subset(mol_pattern_bfp(mol_from_smiles(?1), 2048)) 
                  {'LIMIT ' + str(args.limit) if args.limit is not None else ''}"""
        res = con.execute(sql, (args.query, )).fetchall()

        if args.output is not None:
            sys.stdout = open(args.output, 'wt')
        for smi, i in res:
            print(f'{smi}\t{i}')


if __name__ == '__main__':
    main()
