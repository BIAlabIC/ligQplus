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
from datetime import datetime
import urllib.parse
from requests.adapters import HTTPAdapter 
from requests.packages.urllib3.util.retry import Retry

'''
def search_chembl(chembl_id):
   retry_strategy = Retry(total=3, status_forcelist=[429, 500, 502, 503, 504], 
                       method_whitelist=["HEAD", "GET", "OPTIONS", "POST"])
   adapter = HTTPAdapter(max_retries=retry_strategy)
   http = requests.Session() 
   http.mount("https://", adapter) 
   http.mount("http://", adapter)
   r = http.get(f"https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/{chembl_id}/1")
   if r.ok:
       data =  r.json()
       for x in data:
           x["src_name"] = sources[x["src_id"]]
       return data
   raise Exception(r.text)'''



class Compound_cache():
    headers=['Representative_ligand', 'Zinc_result', 'SMILES', 'Chemmol']
    def __init__(self, cache_path, ligand_dict_path):
        self.cache_path=cache_path
        self.df_cache=None
        self.dict_lig_representative={}
        self.ligand_dict_path=ligand_dict_path
        
    def init(self):
        if os.path.exists(self.cache_path):
            self.df_cache=pd.read_csv(self.cache_path)
            self.df_cache['Chemmol']=[Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)) for smiles in self.df_cache['SMILES']]
        else: 
            #headers=['Representative_ligand', 'Zinc_result', 'SMILES', 'Chemmol']
            self.df_cache=pd.DataFrame(columns=Compound_cache.headers)
        if os.path.exists(self.ligand_dict_path):
            with open(self.ligand_dict_path) as h:
                self.dict_lig_representative=json.load(h)
                

    def find_by_tanimoto(self, smiles, tanimoto):
        self.df_cache['Tanimoto'] = DataStructs.BulkTanimotoSimilarity(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)),
                                                        list(self.df_cache.Chemmol))
        df_cache_filtrado=self.df_cache[self.df_cache['Tanimoto']>=tanimoto]
        if df_cache_filtrado.empty:
            return None
        else:
            df_cache_filtrado=df_cache_filtrado.sort_values('Tanimoto', ascending=False)
            cache_record=df_cache_filtrado.iloc[0]
            return cache_record.to_dict()
        
    def find_by_id(self, ligand_id):
        if ligand_id in self.dict_lig_representative:
            lig_representative=self.dict_lig_representative[ligand_id]
            df_cache_filtrado=self.df_cache[self.df_cache['Representative_ligand']==lig_representative]
            cache_record=df_cache_filtrado.iloc[0]
            return cache_record.to_dict()
        else:
            return None      
        
    
    def put(self, cache_record, cluster_ligands):
        print(cache_record)
        df_cache_record=pd.DataFrame(cache_record)
        self.df_cache=pd.concat([self.df_cache, df_cache_record], ignore_index=True).drop_duplicates('Representative_ligand')
        #print(self.df_cache)
        for cluster_ligand in cluster_ligands:
            self.dict_lig_representative[cluster_ligand]=cache_record['Representative_ligand']
        
    def save(self):
        self.df_cache[['Representative_ligand', 'Zinc_result', 'SMILES']].to_csv(self.cache_path, index=False)
        with open(self.ligand_dict_path, 'w') as h:
            json.dump(self.dict_lig_representative, h)

def get_zinc(smiles, tanimoto_zinc=60):
    '''Busca en zinc ligandos comerciables similares al smiles original con un minimo
    de 60% de similaridad por tanimoto'''
    retry_strategy = Retry(total=3, status_forcelist=[429, 500, 502, 503, 504], 
                        method_whitelist=["HEAD", "GET", "OPTIONS", "POST"])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session() 
    http.mount("https://", adapter) 
    http.mount("http://", adapter)
    smiles_2=urllib.parse.quote(smiles.encode("utf8"))
    request_txt=f"https://zinc20.docking.org/substances.txt?ecfp4_fp-tanimoto-{tanimoto_zinc}={smiles_2}&purchasability=for-sale"
    out=requests.get(request_txt, 
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
        return data
    else:
        sys.stderr.write(f'Fail request to Zinc using SMILES: {request_txt} \n')
        

def purchable_molecules(cluster_dict, dict_ligandid_smiles, cache, tanimoto_cluster, tanimoto_zinc=60): 
    #toma el primer representante de cada cluster para buscar en zinc 
    #print(cache)
    purcheable_per_cluster={}
    for cluster_id, ligand_list in tqdm(cluster_dict.items()): 
        cluster_representative_ligand=ligand_list[0]    
        cluster_representative_smiles=dict_ligandid_smiles[cluster_representative_ligand]
        cache_record=cache.find_by_tanimoto(cluster_representative_smiles, tanimoto_cluster)
        date=datetime.now() #['Representative_ligand', 'Zinc_result', 'SMILES', 'Chemmol']
        # print(cache_record)
        if not cache_record:
            purchable=get_zinc(cluster_representative_smiles, tanimoto_zinc) 
            sleep(0.5)
            if  purchable==None:
                purchable=[]
            dict_cache={'Representative_ligand':cluster_representative_ligand, 
                        'Zinc_result':purchable, 
                        'SMILES':cluster_representative_smiles,
                        'Chemmol':Chem.RDKFingerprint(Chem.MolFromSmiles(cluster_representative_smiles))}
            cache.put(dict_cache, ligand_list)
            purcheable_record={'Date':date.strftime("%d/%m/%Y %H:%M"),
                               'Representative_ligand': cluster_representative_ligand, 
                               'Ligand_list':ligand_list,
                               'Zinc_result':purchable}
        else:
            purcheable_record={'Date':date.strftime("%d/%m/%Y %H:%M"),
                               'Representative_ligand':cluster_representative_ligand, 
                               'Cache_ligand':cache_record['Representative_ligand'],
                               'Ligand_list':ligand_list,
                               'Zinc_result':cache_record['Zinc_result']}
        purcheable_per_cluster[cluster_id]= purcheable_record
    return purcheable_per_cluster

def build_dict_smiles_from_csv(file):
    df_lista_smiles=pd.read_csv(file)
    my_dictionary=dict(zip(df_lista_smiles['ID'],df_lista_smiles['ID_SMILES'].apply(eval)))
    
    return my_dictionary

def fragment_cluster(tupple_ligand_id_smiles_list, dict_ligandid_chemmol, tanimoto_cluster, chunk_size):
    ligs_chunks = [tupple_ligand_id_smiles_list[i:i + chunk_size] 
                   for i in range(0, len(tupple_ligand_id_smiles_list), chunk_size)]
    initial_index=0
    dic_cluster_ligandlist={}
    for ligs_chunk in ligs_chunks:
        dic_cluster_ligandlist_chunk=cluster_maker(ligs_chunk, dict_ligandid_chemmol, tanimoto_cluster)
        for cluster_id, ligand_list in dic_cluster_ligandlist_chunk.items():
            dic_cluster_ligandlist[int(cluster_id)+initial_index]=ligand_list
        initial_index=len(dic_cluster_ligandlist)
    return dic_cluster_ligandlist
        
def cluster_maker(tupple_ligand_id_smiles_list, dict_ligandid_chemmol, tanimoto_cluster): #usar 0.9 en el argparse
    cutoff = 1-tanimoto_cluster
    # print(tupple_ligand_id_smiles_list)
    # print(tanimoto_cluster)
    fps =[]
    for (ligand_id, smiles) in tupple_ligand_id_smiles_list:
        # print(ligand_id, 'so vo')
        # print(dict_ligandid_chemmol[ligand_id], 'soy yo')
        fp=AllChem.GetMorganFingerprintAsBitVect(
                         dict_ligandid_chemmol[ligand_id],2,1024)
        fps.append(fp)
    # fps = [AllChem.GetMorganFingerprintAsBitVect(
    #                 dict_ligandid_chemmol[ligand_id],2,1024)
    #                    for (ligand_id, smiles) in tupple_ligand_id_smiles_list]
    # first generate the distance matrix (Triangular matrix):
    dists = []
    nfps = len(fps)    
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    # now cluster the data:
    #print('-----', dists,nfps,cutoff)
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    dic_cluster_ligandlist=dict()
    for cluster_id, ligand_index_list in enumerate(cs):
        dic_cluster_ligandlist[cluster_id]=[tupple_ligand_id_smiles_list[index][0] for index in ligand_index_list]
    return dic_cluster_ligandlist


def process_lt(lt, tupple_ligand_id_smiles_list, output_path, cache, tanimoto_cluster, chunk_size):
    dict_ligandid_chemmol={}
    dict_ligandid_smiles={}
    tupple_ligand_id_smiles_list_filtered=[]
    #print('---', tupple_ligand_id_smiles_list)
    for (ligand_id, smiles) in tupple_ligand_id_smiles_list:
        chemmol=Chem.MolFromSmiles(smiles)
        if chemmol:
            tupple_ligand_id_smiles_list_filtered.append((ligand_id, smiles))
            assert chemmol, f'No se pudo convertir en chemmol {chemmol}'
            dict_ligandid_chemmol[ligand_id]=chemmol
            dict_ligandid_smiles[ligand_id]=smiles
    #mols=smiles_to_chemmol(lista_smiles_lt)
    cluster_dict=fragment_cluster(tupple_ligand_id_smiles_list_filtered, dict_ligandid_chemmol, tanimoto_cluster, chunk_size)
    zinc_data = purchable_molecules(cluster_dict, dict_ligandid_smiles, cache, tanimoto_cluster) #Busca las moleculas comprables en zinc con hasta tanimoto en 0.6
    return zinc_data

def process_proteome(lt_smiles_dict, output_path, tanimoto_cluster, cache, chunk_size):
    for lt, tupple_ligand_id_smiles_list in tqdm(lt_smiles_dict.items()):
        sys.stderr.write(f'Running lt ID= {lt}\n')
        lt=lt.strip()
        if os.path.isfile(output_path+f'{lt}.csv'): #si esta el archivos quiere decir que ya corrio para ese lt
            sys.stderr.write(f'{lt}.csv already exists\n')
        else:
            if tupple_ligand_id_smiles_list:
                #ver si se puede validar que haya corrido antes un lt similar
                zinc_data=process_lt(lt, tupple_ligand_id_smiles_list, output_path, cache, tanimoto_cluster, chunk_size)
                df_final2=pd.DataFrame(zinc_data)
                df_final=df_final2.transpose()
            else:
                headers=['Cluster_ID', 'Date', 'Representative_ligand', 'Ligand_list', 'Zinc_result']
                df_final=pd.DataFrame(columns=headers)
            df_final.to_csv(path_or_buf=output_path+f'{lt}.csv', index=True)
            
    
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help= 'ORG_lt_id_smile.csv file. Second  \
                        column must be an id-SMILES tupple', required=True)
    parser.add_argument('-o','--output', help= 'Directory path', required=True)
    #parser.add_argument('-dict','--dict', help= 'dictionary uniprot:lt ', required=True)
    parser.add_argument('-p', '--percent', help= 'Identity \
                        between clusters, default=0.9. must be float 0-1', 
                        type=float, default=0.9)
    parser.add_argument('-db', '--database', help= '', default='./cache.csv' ) 
    parser.add_argument('-ld', '--ligand_dict', help= '', default='./cache.json' ) 
    parser.add_argument('-c', '--chunk_size', help= '', 
                        type=int, default=5000)
    args=parser.parse_args()  
    assert os.path.exists(args.input), f'{args.input} does not exist'
    cache=Compound_cache(args.database, args.ligand_dict)
    cache.init()
    try:
        lt_smiles_dict=build_dict_smiles_from_csv(args.input) #.csv lt, [(ligand_id, smiles),...]
        purchable=process_proteome(lt_smiles_dict, args.output, args.percent, cache, args.chunk_size)
        
    except RequestException as e:
        print('Error using API. Try again.')
        raise
    finally:
        cache.save()

if __name__=='__main__':
    main()

