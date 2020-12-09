from evcouplings.utils import read_config_file, write_config_file
from evcouplings.utils.pipeline import execute
import itertools
import os
import argparse
import subprocess

def evcouplingbatch(protid1, protid2, output):
    # opens the basic config file to be modified
    config = read_config_file("/home/cns-mccafferty/Ecoli_CoEv_modeling/basic_config.txt", preserve_order=True)
    
    # specifies output for output pair
    config["global"]["prefix"] = "output/" + output
    
    # sequence id for proteins to be analyzed
    config["align_1"]["sequence_id"] = protid1
    config["align_2"]["sequence_id"] = protid2
    
    # writes out config file
    write_config_file(output + "_config.txt", config)
    print(output + "_config.txt was written.")
    
    # execute config file
    # config = read_config_file(output + "_config.txt")
    # outcfg = execute(**config)
    subprocess.run(["evcouplings_runcfg", output + "_config.txt"])
    

def Main(proteinlist):
    protein_group = list(itertools.combinations(proteinlist,2))


    for i in protein_group:
        output = i[0] + '_' + i[1]
        print('evcoupling for ' + i[0] + ' and ' + i[1] + ' is being analyzed.')
        evcouplingbatch(i[0],i[1],output)
        print('evcoupling for ' + i[0] + ' and ' + i[1] + ' is finished.')


    os.remove(output + "_config.txt")
    print(output + "_config.txt removed.")

parser = argparse.ArgumentParser(description='Calculates pairwaise coevolutionary scores for protein pairs in an interaction group.')
parser.add_argument('proteinlist', help='list of uniprot ids of proteins in a complex groups.', nargs='+')
args = parser.parse_args()

Main(args.proteinlist)
