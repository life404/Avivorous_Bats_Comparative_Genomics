import pandas as pd
import numpy as np
from utils import cmd_run_multiple
from pathlib import Path
import argparse

ftplink = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
def blastdb(genome, output):
    root_path = Path(output)

    db_path = root_path / "genomeIndex"
    db_path.mkdir(exist_ok=True, parents=True)

    mkblastdb_cmd = f"""
	~/TOOLS/ncbi-blast-2.14.0+/bin/makeblastdb -in {genome} -dbtype nucl -out {db_path}/genome_db
    """
    check = cmd_run_multiple(mkblastdb_cmd, output, check="mkblastdb")
    return check, f"{db_path}/genome_db"


def genome_blast(db, output, virus_dir):
    root_path = Path(output)

    output_gag = root_path / "gag"
    output_env = root_path / "env"
    output_pol = root_path / "pol"
    output_gag.mkdir(exist_ok=True, parents=True)
    output_env.mkdir(exist_ok=True, parents=True)
    output_pol.mkdir(exist_ok=True, parents=True)

    tblastn_cmds = [
        f"""
    ~/TOOLS/ncbi-blast-2.14.0+/bin/tblastn -outfmt "6 qseqid sseqid slen sstart send evalue length qseq sseq sframe" -query {virus_dir}/GAG.faa -db {db} -evalue 0.009 -out {output_gag}/genome_blast.outfmt6 -num_threads 64
    """,
        f"""
    ~/TOOLS/ncbi-blast-2.14.0+/bin/tblastn -outfmt "6 qseqid sseqid slen sstart send evalue length qseq sseq sframe" -query {virus_dir}/ENV.faa -db {db} -evalue 0.009 -out {output_env}/genome_blast.outfmt6 -num_threads 64
    """,
        f"""
    ~/TOOLS/ncbi-blast-2.14.0+/bin/tblastn -outfmt "6 qseqid sseqid slen sstart send evalue length qseq sseq sframe" -query {virus_dir}/POL.faa -db {db} -evalue 0.009 -out {output_pol}/genome_blast.outfmt6 -num_threads 64
    """,
    ]

    check = cmd_run_multiple(tblastn_cmds, ou_path=root_path, check="genome_tblastn")

    return check, [output_gag, output_env, output_pol]


def reverse_blast(output, virus_db, results_dir):
    root_path = Path(output)

    tblastn_cmds = list()
    results_path = list()
    for dir in results_dir:
        if dir.name in ['env', 'gag']:
            len_threshold = 200
        elif dir.name == 'pol':
            len_threshold = 400
        else:
            print('something wrongs')
            exit(1)
            
        tblastn_result = Path(dir) / "genome_blast.outfmt6"
        output = Path(dir) / "virus_blast.outfmt6"
        results_path.append(output)
        first_blast_seq = parse_genome_tblastn(tblastn_result, len_threshold)
        tblastn_cmds.append(
            f"""~/TOOLS/ncbi-blast-2.14.0+/bin/tblastn -outfmt "6 qseqid sseqid slen qstart qend evalue length qseq sseq sframe" -query {first_blast_seq} -db {virus_db} -evalue 1e-5 -out {output} -num_threads 64"""
        )
    check = cmd_run_multiple(tblastn_cmds, ou_path = root_path, check = 'virus_tblastn')
    return results_path


def parse_genome_tblastn(tblastn_result, len_threshold):
    results = pd.read_csv(
        tblastn_result,
        sep="\t",
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "slen",
            "sstart",
            "send",
            "evalue",
            "length",
            "qseq",
            "sseq",
            "sframe",
        ],
    )
    results = results[results.length > len_threshold]
    first_blast_seq_path = tblastn_result.parent / "first_blast_hit_seq.fa"
    blasted_probs = open(first_blast_seq_path, "w")
    for index, result in results.iterrows():
        if result.sstart < result.send:
            direction = "P"
        else:
            direction = "N"
        blasted_probs.write(
            f""">{result.sseqid}-{result.sstart}-{result.send}-{result.length}-{direction}\n{result.sseq.strip().replace('-','')}\n"""
            )
    blasted_probs.close()
    return first_blast_seq_path

def parse_virus_tblastn(tblastn_result):
    results = pd.read_csv(
        tblastn_result,
        sep='\t',
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "slen",
            "qstart",
            "qend",
            "evalue",
            "length",
            "qseq",
            "sseq",
            "sframe", 
        ]
    )
    tmp_df = results.qseqid.str.split('-', expand=True)
    tmp_df.columns = ['scaffold', 'firststart', 'firstend', 'firstlen', 'firstdirection']
    results = pd.concat([tmp_df, results], axis=1)
    results = results.drop_duplicates('qseqid', keep='first')
    results.firststart = results.firststart.apply(int)
    results.firstlen = results.firstlen.apply(int)
    results.firstend = results.firstend.apply(int)
    results = results[results.length > results.firstlen * 0.7]
    return results

def merge_results(df):
    if len(df) == 1:
        return df.loc[:, ['firststart', 'firstend']]
    else:
        df = df.sort_values('firststart').reset_index(drop=True)
        hit = df.loc[0, :]
        hit_range = np.arange(min(hit.firststart, hit.firstend), max(hit.firststart, hit.firstend)+1)
        hit_ranges = {'firststart': [],
                      'firstend': [],
                    } 
        for i in range(1, len(df)):
            tmp = df.loc[i, :]
            tmp_range = np.arange(min(tmp.firststart, tmp.firstend), max(tmp.firststart, tmp.firstend)+1)
            if np.intersect1d(hit_range, tmp_range).size > 0:
                hit_range = np.union1d(hit_range, tmp_range)
                hit_ranges['firststart'].append(min(hit_range))
                hit_ranges['firstend'].append(max(hit_range))
            else:
                if 0 < max(hit_range) - min(tmp_range) <= 50:
                    hit_range = np.arange(min(hit_range), max(tmp_range)+1)
                    hit_ranges['firststart'].append(min(hit_range))
                    hit_ranges['firstend'].append(max(hit_range))
                else:
                    hit_ranges['firststart'].append(min(hit_range))
                    hit_ranges['firstend'].append(max(hit_range))
                    hit_ranges['firststart'].append(min(tmp_range))
                    hit_ranges['firstend'].append(max(tmp_range))
                    hit_range = tmp_range
        final_df = pd.DataFrame(hit_ranges).sort_values('firststart').drop_duplicates('firststart', keep='first').drop_duplicates('firstend', keep='last')
        return pd.DataFrame(final_df).reset_index(drop=True)

def final_probs(results):
    return results.groupby(['sseqid', 'scaffold', 'firstdirection']).apply(merge_results).reset_index()
    
def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-g', '--genome', help='The path of genome file')
    parse.add_argument('-v', '--virus_dir', help='The directory contains threee [gag, env, pol] retrovirus protein fasta files')
    parse.add_argument('--virus_db', default='/home/panda2bat/.ervin/virus_db_store/Viruses', help=f'The ncbi nuclitaed database of lastest Viruses squences from NCBI, the ftp link is:\n{ftplink}') 
    parse.add_argument('-o', "--output", default='The output directory')
    args = parse.parse_args()
    return args 

def main():
    args = make_parse()
    genome = args.genome 
    output = args.output
    virus_dir = args.virus_dir
    virus_db = args.virus_db
    virus_info = pd.read_csv("/home/panda2bat/.ervin/virus_db_store/virus.info.tsv", header=None, names=['sseqid', 'species'], sep='\t')

    check, db = blastdb(genome, output)
    check, results_dir = genome_blast(db, output, virus_dir)
    reverse_blast_results = reverse_blast(output, virus_db=virus_db, results_dir=results_dir)
    for result_file in reverse_blast_results:
        results_dir = result_file.parent
        result = parse_virus_tblastn(result_file)
        final_tbl = final_probs(result)
        final_tbl = final_tbl.drop(columns='level_3')
        final_tbl.to_csv((results_dir / f"{result_file.parts[-2]}.tsv"), sep='\t', index=False)
        pd.merge(final_tbl, virus_info, on='sseqid', how='inner').groupby('sseqid').size().to_csv((results_dir / f"{result_file.parts[-2]}.summary"), header=['Counts'], sep='\t')
        
        
if __name__ == "__main__":
    main()
