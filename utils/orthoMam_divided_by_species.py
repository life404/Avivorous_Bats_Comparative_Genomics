#! /usr/bin/python3

import argparse
from pathlib import Path
from Bio import SeqIO

parse = argparse.ArgumentParser()
parse.add_argument('--input', dest='input', type = str, help = 'The input directory of input')
parse.add_argument('--output', dest = 'output', type = str, help = 'The output directory of results')
args = parse.parse_args()

input_d = args.input
output = args.output

fa_lst = list(Path(input_d).glob('*.fasta'))

species_lst = {}
for fa in fa_lst:
    records = SeqIO.parse(fa, 'fasta')
    for record in records:
        species = record.id
        if species in species_lst.keys():
            record.id = fa.stem.replace('_AA','').replace('_NT', '')
            species_lst[species].append(record)
        else:
            species_lst[species] = []
            record.id = fa.stem.replace('_AA','').replace('_NT', '')
            species_lst[species].append(record)

Path(output).mkdir(exist_ok=True, parents=True)
for specie, records in species_lst.items():
    SeqIO.write(records, f'{output}/{specie}.fa', 'fasta')
    
    


        