import argparse
from collections import defaultdict


def PFAM_ID(input_file):
    """file with pfam accession list"""
    pfam_ids = open(input_file)
    lines=pfam_ids.readlines()
    strip_pfam_ids = [line.strip().split('.')[0] for line in lines ]

    pfam_ids.close()
    return strip_pfam_ids

def pdb_pfam_mapping(input_file):
    """ pdb_pfam_mapping
    columns: 0=PDB_ID; 1=CHAIN_ID; 2=RES_NUM_START; 3=RES_NUM_STOP; 4=PFAM_ACC """
    mapping = open(input_file)
    mapping_lines=mapping.readlines()
    mapping.close()
    dic_out= defaultdict(list)
    for line in mapping_lines:
        split_mapping_lines=line.split("\t")
        pfam = split_mapping_lines[4].split('.')[0]
        dominio_pdb={'pdb':split_mapping_lines[0], 'chain':split_mapping_lines[1],
        'start':split_mapping_lines[2], 'stop':split_mapping_lines[3]}
        dic_out[pfam].append(dominio_pdb)


    return dic_out

def id_cross(PFAM_data, mapping_data):
    """Cross IDs between inputs file"""

    for pfam in PFAM_data:
        if pfam in mapping_data:
            for pdb in mapping_data[pfam]:
                print(pfam+"\t"+ f'{pdb["pdb"]}\t{pdb["chain"]}\t{pdb["start"]}\t{pdb["stop"]}')
    return 0


def parse_arguments():

    parser = argparse.ArgumentParser(description='Extract PDB from PFAM IDs')
    parser.add_argument("-l", '--PFAM_list', default=None, help="Input file should have PFAM ID list")
    parser.add_argument("-m", '--mapping_pdb_pfam',  default="pdb_pfam_mapping.txt", help="Input file of PDB and PFAM ID link: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt ")

    return parser

def main():

    parser=parse_arguments()
    args=parser.parse_args()

    id_cross(PFAM_ID(args.PFAM_list), pdb_pfam_mapping(args.mapping_pdb_pfam))
    return 0


if __name__=='__main__':
    main()
