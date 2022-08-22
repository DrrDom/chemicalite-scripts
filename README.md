### The repository of useful scripts to manipulate SQLite DB using Chemicalite

Normally you should create DB from SMILES, add the mol field and add desired fingerprints:  
`create_db.py` - create DB from SMILES file  
`create_mol_field.py` - create Mol objects in DB  
`add_fp_to_db.py` - add chosen fingerprints to DB - morgan, atom-pairs or pattern (for substructure search) are supported.

After that you may perform similarity of substructure search:  
`sim_search.py` - perform similarity search using the chosen fingerprints  
`substr_search.py` - perform substructure search (after pattern fp were added)  

Examples:  
```
create_db.py -i input.smi -o database.db
create_mol_field.py -d database.db
add_fp_to_db.py -d database.db -f pattern
add_fp_to_db.py -d database.db -f morgan

sim_search.py -d database.db -q 'CCO' -f morgan
substr_search.py -d database.db -q 'CCO' 
```

##### Dependency

`rdkit` - https://www.rdkit.org/  
`Chemicalite` - SQLite cartridge for RDKit - https://github.com/rvianello/chemicalite
