#!/usr/bin/env python3

import argparse
import sqlite3


def main():
    parser = argparse.ArgumentParser(description='Similarity search in SQLite DB.')
    parser.add_argument('-d', '--input_db', metavar='FILENAME', required=True,
                        help='input SQLite DB.')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where Mol objects are stored. Default: mols.')
    parser.add_argument('-m', '--mol_field', metavar='STRING', default='mol',
                        help='column which contains Mol object. Default: mol.')
    parser.add_argument('-f', '--fp', metavar='STRING', default='morgan', choices=['morgan', 'pattern', 'atom_pairs'],
                        help='fingerprint type to compute. Default: morgan.')

    args = parser.parse_args()
    with sqlite3.connect(args.input_db) as con:

        con.enable_load_extension(True)
        con.load_extension('chemicalite')
        con.enable_load_extension(False)

        # create a virtual table to be filled with bfp data
        # con.execute(f"DROP TABLE IF EXISTS {args.table}_{args.fp}_idx")
        con.execute(f"CREATE VIRTUAL TABLE IF NOT EXISTS {args.table}_{args.fp}_idx "
                    f"USING rdtree(id, fp bits(2048))")

        # compute and insert the fingerprints
        sql = f"INSERT INTO {args.table}_{args.fp}_idx(id, fp) " \
              f"SELECT rowid, mol_{args.fp}_bfp({args.mol_field}, {'2,' if args.fp == 'morgan' else ''} 2048) " \
              f"FROM {args.table} " \
              f"WHERE {args.mol_field} IS NOT NULL"
        con.execute(sql)

        con.commit()


if __name__ == '__main__':
    main()
