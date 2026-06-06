#! /usr/bin/python3

import argparse
from pathlib import Path
import shutil

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', dest = 'input', help = 'The input directory or a single input file, if it is a directory, it should only contain fasta files')
    parse.add_argument('-o', dest = 'output', help = 'If input is a directory, the paragram will generate a result directory for each input. If input is a file, the output parameter will be used as output directory for this input')
    parse.add_argument('-t', dest = 'tree', help = 'The labeled tree for codeml')
    parse.add_argument('-c', dest = 'config', help = 'The raw config file')
    parse.add_argument('-s', dest = 'species', default = 'null', help = 'The front ground species used to label the tree, if the number of species is one, just input the species, if two, input just like (A, B)')
    args = parse.parse_args()
    return(args)

def codeml_commands_single(config, output, tree, seq, species):
    root = Path(output)
    root.mkdir(exist_ok=True, parents=True)

    tree_path = output / Path(tree).name
    tree_content = open(tree, 'r').read()
    tree_condtent = tree_content.replace(species, f"{species} #1")
    tree_tmp = open(tree_path, 'w')
    tree_tmp.writelines(tree_condtent)
    tree_tmp.close()
    
    #shutil.copy(tree, tree_path)

    seq_path = output / Path(seq).name
    #seq_path= f"./{Path(seq).name}"
    shutil.copy(seq, seq_path)

    config_path = root / 'codeml.ctl'
    config_content = open(config, 'r').read()
    config_content = config_content.replace('SEQZ', f"./{seq_path.name}")
    config_content = config_content.replace('TREEZ', f"./{tree_path.name}")
    config_tmp = open(config_path, 'w')
    config_tmp.writelines(config_content)
    config_tmp.close()

    cmd = f"cd {root};codeml codeml.ctl"
    #check = output / 'codeml.ok'
    #return_code, output_info, command = cmd_run(cmd, ou_path = root)
    #if return_code == 0:
    #    check.touch()
    #else:
    #    print(output_info)
    #    exit(1)
    print(cmd)

#def codeml_commands_multiple

def main():
    args = make_parse()
    input = Path(args.input)
    output = Path(args.output)
    tree = Path(args.tree)
    config = Path(args.config)
    species = args.species
    
    if input.is_file():
        seq = input 
        codeml_commands_single(config, output, tree, seq, species)
    elif input.is_dir():
        pass
        

if __name__ == '__main__':
    main()