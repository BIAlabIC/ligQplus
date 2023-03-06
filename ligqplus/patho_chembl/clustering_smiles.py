#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 10:20:51 2021

@author: fleer
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
from rdkit.ML.Cluster import Butina
from collections import defaultdict
import pandas as pd
import requests
import argparse
import os
import sys
import json
from requests.exceptions import RequestException 
from tqdm import tqdm
from time import sleep
from contextlib import closing

def ClusterFps(fps, similaridad):
    cutoff = 1-similaridad
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

def process_csv(file):
    df_lista_smiles=pd.read_csv(file)
    return df_lista_smiles

def csv_to_dict(file):
    df_lista_smiles= process_csv(file)
    my_dictionary=dict(zip(df_lista_smiles['ID'],df_lista_smiles['ID_SMILES'].apply(eval)))
    return my_dictionary


def smiles_to_chemmol(lista_smiles_lt):    
    '''transform csv file into a list of [[chem.molecule, molecule_id]] using RDKit'''
    mols = [[Chem.MolFromSmiles(mol[1]),mol[0]] for mol in lista_smiles_lt ]
    #print(mols)
    return mols

#cluster por SMILES
def save_as_smiles(clusters, mols, dic):
    clusters_smiles = []
    for cluster in clusters:
        smiles = []
        for molec in cluster:
            for lig in dic:
                smile_id=mols[molec][1]
                if lig[0]==smile_id:
                       smiles.append(lig[1])
        if len(smiles)>=1:    
            clusters_smiles.append(smiles)
    clusters_smiles.sort(key=len,reverse = True)
    return clusters_smiles

def save_as_id(clusters, mols):
    clusters_smiles = []
    for cluster in clusters:
        smiles = []
        for molec in cluster:
            smiles.append(mols[molec][1])
        if len(smiles)>1:
            clusters_smiles.append(smiles)
    clusters_smiles.sort(key=len,reverse = True)
    return clusters_smiles

def save_as_tuple(clusters, mols, dic):
    clusters_tuple = []
    for cluster in clusters:
        smiles = []
        for molec in cluster:
            for lig in dic:
                smile_id=mols[molec][1]
                if lig[0]==smile_id:
                       smiles.append((smile_id, lig[1]))
        if len(smiles)>1:    
            clusters_tuple.append(smiles)
    clusters_tuple.sort(key=len,reverse = True)
    return clusters_tuple



def cluster_maker(mols, tuple_list, percent, save_as='smiles'): #usar 0.9 en el argparse, save_as='smiles'
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol[0],2,1024) for mol in mols if mol[0] != None]
    clusters = ClusterFps(fps, percent)
    smiles_lt=tuple_list
    if save_as=='smiles':
        smiles_list=save_as_smiles(clusters, mols, smiles_lt)
        return smiles_list
    elif save_as=='id_list':
        id_list=save_as_id(smiles_lt, mols)
        return id_list
    elif save_as=='both':
        tuple_id_smiles=save_as_tuple(clusters, mols, smiles_lt)
        return tuple_id_smiles
    


def get_zinc(smiles, tanimoto=60):
    text=[]
    out=requests.get(f"http://zinc20.docking.org/substances.txt?ecfp4_fp-tanimoto-{tanimoto}={smiles}&purchasability=for-sale", 
                    stream=True, allow_redirects=True)
    if out.status_code==200:
        with closing(out) as r, open('/tmp/zinc', 'w') as h:
            for x in r.iter_content(chunk_size=512*1024):
                h.write(x.decode())
        with open('/tmp/zinc', 'r') as h:
            info=h.readlines()
        data=[]
        if info:
            if len(info)>1 and info!=None:
                for zinc_smiles in info:
                    zinc_id=zinc_smiles.split('\t')[0]
                    smiles_id=zinc_smiles.split('\t')[1]
                    if zinc_id:
                        data.append((zinc_id, smiles_id))
                #print(data)
        #else:
         #   pass
            #zinc_id=[]
            #data.append(zinc_id)
        return data

def purchable_molecules(clusters, tanimoto=60): 
    pur_per_cluster=[]
    for cluster in clusters:
        purch_smiles=[]
        for smiles in cluster:
            purchable=get_zinc(smiles, 60)
            sleep(0.5)
            purch_smiles.append(purchable)
        pur_per_cluster.append(purch_smiles)
    return pur_per_cluster


def process_lt_dict(lt_file):
    with open(lt_file) as h: 
        lt_dict = json.load(h)
    return lt_dict


def search_per_cluster(file, path, lt_uniprot_dict, percent=0.9): #TODO
    my_dictionary=csv_to_dict(file)
    lista_smiles_lt = []
    lt_dict=process_lt_dict(lt_uniprot_dict)
    for k, v in tqdm(lt_dict.items()):
        lt_dict={}
        sys.stderr.write(f'Running lt ID= {v}\n')
        filepath = path
        if os.path.isfile(filepath+f'{v}.csv'):
            print(f'{v}.csv already exists\n')
            pass
        else:
            for lt, ligans in my_dictionary.items():
                lt=lt.strip()
                data_lt=[]
                if lt == v:
                    if ligans!=[]:
                        for item in ligans:
                            lista_smiles_lt.append(item)
                        mols=smiles_to_chemmol(lista_smiles_lt)
                        cluster=cluster_maker(mols, lista_smiles_lt, percent)
                        zinc_data = purchable_molecules(cluster)
                        #print(zinc_data)
                        for cluster_id, cluster_list in enumerate(zinc_data): 
                            for ligand_zinc in cluster_list:
                                #print(ligand_zinc)
                                if ligand_zinc:
                                    for tup in ligand_zinc:
                                        if tup[0]:
                                            zinc_id=tup[0]
                                            zinc_smiles=tup[1]
                                            data_lt.append({'cluster_id':cluster_id, \
                                                                      'zinc_id':zinc_id, 'zinc_smiles':zinc_smiles})
                                else:
                                    data_lt.append({'cluster_id':cluster_id, \
                                                                      'zinc_id':"no zinc ID", 'zinc_smiles':"-"})
                    else:
                        print(f'{v} has no zinc ids')
                    df_final=pd.DataFrame(data_lt)
                    df_final.to_csv(path_or_buf=f'{filepath}/{v}.csv', index=False)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help= 'ORG_lt_id_smile.csv file. Second  \
                        column must be an id-SMILES tupple', required=True)
    parser.add_argument('-o','--output', help= 'Directory path', required=True)
    parser.add_argument('-dict','--dict', help= 'dictionary uniprot:lt ', required=True)
    parser.add_argument('-p', '--percent', help= 'Identity \
                        between clusters, default=0.9. must be float 0-1', default=0.9)
    args=parser.parse_args()
    try:
        purchable=search_per_cluster(args.input, args.output, args.dict, args.percent)
    except RequestException as e:
        print('Error using API. Try again.')
        raise

if __name__=='__main__':
    main()
