#!/usr/bin/env python3

import argparse
import sqlite3
import sys


def main():
    parser = argparse.ArgumentParser(description='Create a DB with a table containing id and smi columns.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, default=None,
                        help='input SMILES file. The file should contain names which should be distinct, records with '
                             'repeated ids will be ignored during insert.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='SQLite DB to create.')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where to store smiles and ids of molecules. Default: mols.')
    parser.add_argument('-s', '--sep', metavar='STRING', default=None,
                        help='separator in input SMILES file. Default: None (whitespaces).')
    parser.add_argument('-v', '--verbose', required=False, default=False, action='store_true',
                        help='print progress to stderr.')

    args = parser.parse_args()
    with sqlite3.connect(args.output) as con:

        con.execute(f"CREATE TABLE {args.table} (id TEXT UNIQUE, smi TEXT)")
        con.commit()

        i = 0
        lines = []
        with open(args.input) as f:
            for line in f:
                lines.append(line.strip().split(args.sep))
                i += 1
                if i % 10000 == 0:
                    con.executemany(f"INSERT OR IGNORE INTO {args.table} (smi, id) VALUES(?, ?)", lines)
                    con.commit()
                    lines = []
                    if args.verbose:
                        sys.stderr.write(f'\r{i} records were processed')
            con.executemany(f"INSERT OR IGNORE INTO {args.table} (smi, id) VALUES(?, ?)", lines)
            con.commit()
            if args.verbose:
                sys.stderr.write(f'\r{i} records were processed')


if __name__ == '__main__':
    main()
