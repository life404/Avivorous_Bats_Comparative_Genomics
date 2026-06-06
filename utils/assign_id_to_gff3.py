#! /usr/bin/python3

import pandas as pd
import numpy as np
import re
import argparse
from pathlib import Path


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "--gff3",
        default="/home/panda2bat/Avivorous_bat/output/08_genome-annotation/N.aviator/PASA_update/test.gff3",
        help="The path of gff3",
    )
    parse.add_argument(
        "--uniprot",
        default="/home/panda2bat/Avivorous_bat/output/08_genome-annotation/N.aviator/Function/UniProt/blastp.uniprot.out",
        help="The result of blastp against to UniProt",
    )
    parse.add_argument(
        "--uniprot2id",
        default="/mnt/Extra_storage/UniProt/uniprotID_to_geneID",
        help='The file has two columns and delimeted by tab, the first is the UniProt ID, the second column is genename in UniProt, for example "Q6GZX0  005R_FRG3G"',
    )
    args = parse.parse_args()
    return args


def parse_uniprot(uniprot_file, uniprot2id_file):
    uniprot = pd.read_csv(uniprot_file, header=None, index_col=None, sep="\t")
    uniprot = uniprot.loc[:, [0, 1]]
    uniprot.columns = ["evm", "database"]
    uniprot2id = pd.read_csv(uniprot2id_file, header=None, index_col=None, sep="\t")
    uniprot2id.columns = ["database", "gene"]
    uniprot = pd.merge(uniprot, uniprot2id, how="left", on="database")
    uniprot = uniprot.set_index('evm')
    uniprot = uniprot.to_dict()
    return uniprot


def update_gff(gff3, uniprot, nr):
    geneid = 0
    with open(gff3, "r") as file:
        for line in file:
            line = line.strip()
            array = line.split("\t")
            if array[2] == "gene":
                geneid += 1
                evmid = array[8].split(';')[0].replace('ID=', '')
                if evmid in uniprot['database'].keys():
                    line = line + f";Pseudo_id=NAVI{geneid:011d};Mapping_id={uniprot['database'][evmid]};Gene={uniprot['gene'][evmid]}"
                elif evmid in nr['database'].keys():
                    line = line = f";Pseudo_id=NAVI{geneid:011d};Mapping_id={nr['database'][evmid]};Gene={nr['gene'][evmid]}"
                else:
                    line = line = f";Pseudo_id=NAVI{geneid:011d};Mapping_id={nr['database'][evmid]};Gene={nr['gene'][evmid]}" 

                mrna, exon, cds, utr5, utr3 = 1
            elif array[2] == "mRNA":
                line = re.sub(r"Parent=.*", f"Parent=NAVI{geneid:011d}", line)
                line = re.sub(r"ID=.*;", f"ID=NAVI{geneid:011d}.mrna{mrna};", line)
                mrna += 1
            elif array[2] == "CDS":
                line = re.sub(r"Parent=.*", f"Parent=NAVI{geneid:011d}", line)
                line = re.sub(r"ID=.*;", f"ID=NAVI{geneid:011d}.cds{mrna};", line)
                cds += 1
            elif array[2] == "exon":
                line = re.sub(r"Parent=.*", f"Parent=NAVI{geneid:011d}", line)
                line = re.sub(r"ID=.*;", f"ID=NAVI{geneid:011d}.exon{mrna};", line)
                exon += 1
            elif array[2] == "five_prime_UTR":
                line = re.sub(r"Parent=.*", f"Parent=NAVI{geneid:011d}", line)
                line = re.sub(r"ID=.*;", f"ID=NAVI{geneid:011d}.utr5{mrna};", line)
                utr5 += 1
            elif array[2] == "three_prime_UTR":
                line = re.sub(r"Parent=.*", f"Parent=NAVI{geneid:011d}", line)
                line = re.sub(r"ID=.*;", f"ID=NAVI{geneid:011d}.utr3{mrna};", line)
                utr3 += 1
            else:
                continue
            print(line)
