#!/usr/bin/env python3

import argparse
import sqlite3


def main():
    parser = argparse.ArgumentParser(description='Insert a column with RDKit Mol object.')
    parser.add_argument('-d', '--input_db', metavar='FILENAME', required=True,
                        help='input SQLite DB.')
    parser.add_argument('-t', '--table', metavar='STRING', default='mols',
                        help='table name where source data are stored and where generated Mol objects will stored. '
                             'Default: mols.')
    parser.add_argument('-i', '--source_field', metavar='STRING', default='smi',
                        help='source field which will be converted to Mol object. Default: smi.')
    parser.add_argument('-f', '--source_field_type', metavar='STRING', default='smi', choices=['smi', 'molblock'],
                        help='type of input field data. Can be smi for SMILES and mol for MolBlock. Default: smi.')
    parser.add_argument('-o', '--output_field', metavar='STRING', default='mol',
                        help='output field where RDKit Mol objects will be stored. Default: mol.')

    args = parser.parse_args()
    with sqlite3.connect(args.input_db) as con:

        con.enable_load_extension(True)
        con.load_extension('chemicalite')
        con.enable_load_extension(False)

        columns = list(i[1] for i in con.execute(f"PRAGMA table_info({args.table})"))
        if args.output_field not in columns:
            con.execute(f"ALTER TABLE {args.table} ADD COLUMN {args.output_field} MOL")
        if args.source_field_type == 'smi':
            con.execute(f"UPDATE {args.table} "
                        f"SET {args.output_field} = mol_from_smiles({args.source_field}) "
                        f"WHERE {args.output_field} IS NULL")
        elif args.source_field_type == 'molblock':
            con.execute(f"UPDATE {args.table} "
                        f"SET {args.output_field} = mol_from_molblock({args.source_field}) "
                        f"WHERE {args.output_field} IS NULL")
        con.commit()


if __name__ == '__main__':
    main()
