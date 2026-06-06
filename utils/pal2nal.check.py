#! /usr/bin/python3

# This script is used to check the cDNA fasta file and corresponding protein file before runing pal2nal
import sys

sys.path.append("/home/panda2bat/Avivorous_bat/script/")

import argparse
from Bio import SeqIO
from pathlib import Path
import shutil
import concurrent.futures as cf
from utils import cmd_run_multiple


def write_tmp(key, tmp):
    filename = key.split("|")[0]
    SeqIO.write(fna_records[key], f"{tmp}/{filename}.fna", "fasta")
    SeqIO.write(faa_records[key], f"{tmp}/{filename}.faa", "fasta")
    check_cmd = f"~/TOOLS/pal2nal.v14/pal2nal.pl {tmp}/{filename}.faa {tmp}/{filename}.fna > {tmp}/{filename}.pal2nal"
    return check_cmd


parse = argparse.ArgumentParser()
parse.add_argument("-c", "--cdna", type=str, help="The cdna file")
parse.add_argument("-p", "--protein", type=str, help="The protein file")
parse.add_argument(
    "-t", "--tmp", type=str, help="The tmp directory", default=(Path().cwd() / ".tmp")
)
parse.add_argument(
    "-r",
    "--remove",
    type=bool,
    default=False,
    help="if a pair of cdna and protein seuqneces can pass the check of pal2nal, this pair of sequences will be removed from the final results",
)
args = parse.parse_args()

fna_records = {record.id: record for record in SeqIO.parse(args.cdna, "fasta")}
faa_records = {record.id: record for record in SeqIO.parse(args.protein, "fasta")}

sort_fna = sorted(list(fna_records.keys()))
sort_faa = sorted(list(faa_records.keys()))
filename_suffix = list(fna_records.keys())[0].split("|")[1]

if sort_faa != sort_fna:
    print("The id in {args.cdna} is different with {args.protein}")
    exit(1)
else:
    tmp = Path(args.tmp)
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(exist_ok=True, parents=True)

    cmd_lst = list()
    with cf.ThreadPoolExecutor(max_workers=64) as e:
        process_lst = [e.submit(write_tmp, key, tmp) for key in fna_records.keys()]
        for process in cf.as_completed(process_lst):
            cmd = process.result()
            cmd_lst.append(cmd)

    check = cmd_run_multiple(cmd_lst, ou_path=tmp, check="pal2nal.check", num=128)

    if check.exists():
        pal2nal_lst = tmp.glob("*.pal2nal")
        zero_lst = [i for i in pal2nal_lst if i.stat().st_size == 0]
        if len(zero_lst) == 0:
            print(f"All sequences pass the check of PAL2NAL in {filename_suffix}")
        else:
            print(f"Following sequences of {filename_suffix} can not pass the check of PAL2NAL:\n")
            for i in zero_lst:
                print(i)
            if args.remove:
                zero_key = [
                    i.name.replace(".pal2nal", "|" + filename_suffix) for i in zero_lst
                ]
                filtered_faa = [
                    faa_records[key]
                    for key in faa_records.keys()
                    if not key in zero_key
                ]
                filtered_fna = [
                    fna_records[key]
                    for key in fna_records.keys()
                    if not key in zero_key
                ]
                SeqIO.write(
                    filtered_faa, f"{Path(args.protein).absolute()}.filtered", "fasta"
                )
                SeqIO.write(
                    filtered_fna, f"{Path(args.cdna).absolute()}.filtered", "fasta"
                )
                print(
                    f"These {len(zero_key)} sequences have beed remove from final filtered file"
                )
