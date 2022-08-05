#!/usr/bin/env python3

import argparse
import sqlite3


def compose_index_table_name(main_table_name, fp, radius_morgan):
    if fp == 'morgan':
        index_table_name = f'{main_table_name}_{fp}{radius_morgan}_idx'
    else:
        index_table_name = f'{main_table_name}_{fp}_idx'
    return index_table_name


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
    parser.add_argument('-r', '--radius_morgan', metavar='INTEGER', default=2, type=int,
                        help='radius of Morgan fingerprint. Default: 2.')

    args = parser.parse_args()
    with sqlite3.connect(args.input_db) as con:

        con.enable_load_extension(True)
        con.load_extension('chemicalite')
        con.enable_load_extension(False)

        # create a virtual table to be filled with bfp data
        # con.execute(f"DROP TABLE IF EXISTS {args.table}_{args.fp}_idx")
        index_table_name = compose_index_table_name(args.table, args.fp, args.radius_morgan)
        con.execute(f"CREATE VIRTUAL TABLE IF NOT EXISTS {index_table_name} "
                    f"USING rdtree(id, fp bits(2048))")

        # compute and insert the fingerprints
        sql = f"INSERT INTO {index_table_name}(id, fp) " \
              f"SELECT rowid, mol_{args.fp}_bfp({args.mol_field}, {f'{args.radius_morgan},' if args.fp == 'morgan' else ''} 2048) " \
              f"FROM {args.table} " \
              f"WHERE {args.mol_field} IS NOT NULL"
        con.execute(sql)

        con.commit()


if __name__ == '__main__':
    main()
