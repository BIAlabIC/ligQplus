#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates csv file with ligands SMILES from each locus tag in an organism. 
"""
import os
import argparse
import pandas as pd
from tqdm import tqdm
import sys
from ast import literal_eval
import json
import csv

def open_directory(path):
    dire=[]
    base = path
    for subdir in os.listdir(base):
        dire.append(subdir)
    return dire

def process_lt_dict(lt_file):
    with open(lt_file) as h: 
        lt_dict = json.load(h)
    return lt_dict

def process_smiles_dict(smiles_path):
    with open(smiles_path) as dictionary:
        data=dictionary.read()
        smiles_dict = literal_eval(data)
        return smiles_dict
    

def final_csv_maker(path, lt_file, name='ORG_lt_smiles.csv'):
    lt_dict=process_lt_dict(lt_file)
    dire=open_directory(path)
    final_list=[]
    field_names = ['ID', 'SMILES']
    for k, v in tqdm(lt_dict.items()):
        lt_dict={}
        sys.stderr.write(f'Running Uniprot ID= {k}\n')
        filepath = path + os.sep + k
        if not os.path.isfile(filepath+ os.sep +'smiles_dict_filtered'):
            sys.stderr.write(f'{k} does not have smiles_dict_filtered\n')
            continue
        else:
            smiles_dict=process_smiles_dict(filepath+ os.sep +'smiles_dict_filtered')
            smiles_list=[]
            for ids, smiles in smiles_dict.items():
                smiles_list.append(smiles)
            lt_dict['ID']=v
            lt_dict['SMILES']=smiles_list
        final_list.append(lt_dict)
    df_final_list=pd.DataFrame(final_list)
    df_final_list=df_final_list.groupby('ID')
    df_final_list2=df_final_list['SMILES'].sum().apply(set).apply(list)
    df_final_list2=df_final_list2.reset_index()
    df_final_list2.to_csv(path_or_buf=f'{path}/{name}', index=False)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help= 'Directory path', required=True)
    parser.add_argument('-dict', '--dataset', help='Path to the locus tag-uniprot dictionary', required=True)
    parser.add_argument('-o', '--output', help= 'type: ORG_lt_smiles.csv Change ORG.', required=True)
    args=parser.parse_args()
    
    org_lt_ligand=final_csv_maker(args.input, args.dataset, args.output)



if __name__=='__main__':
    main()
