#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 09:21:40 2021

@author: fleer
"""

import pandas as pd
import requests
import os
from sqlalchemy import create_engine
import pandas as pd
import sqlite3
from io import StringIO
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdchem
from rdkit import RDConfig
import argparse
import sys
from tqdm import tqdm
from functools import reduce
from requests.adapters import HTTPAdapter 
from requests.packages.urllib3.util.retry import Retry
#path='/home/fleer/Desktop/Fede/analisis_result'

def open_directory(path):
    dire=[]
    base = path
    for subdir in os.listdir(base):
        dire.append(subdir)
    return dire 

def cleanlist(lis):
    cleanedList = [x for x in lis if str(x) != 'nan']
    return cleanedList

def pdb_to_smiles(ligands, mock=False):
    new_data = []
    if mock:
        return new_data
    retry_strategy = Retry(total=3, status_forcelist=[429, 500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session() 
    http.mount("https://", adapter) 
    http.mount("http://", adapter)  
    r = http.post("https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/",data=",".join(ligands))
   # print(r)
    if r.ok:
        data =  r.json()  
       # print(data)

        for k,v in data.items():
            r = v[0]            
            r2 = {"pdb_ligand":k}
            if len(r["smiles"]):
                r2["smiles"] = r["smiles"][0]["name"]           
                new_data.append(r2) 
        return new_data
    else:
        print(r)
        return new_data
   # raise Exception(r.text)
    
def pdb_smiles_dic(pdb_list, n=50, mock=False):
    pdb_smiles=[]
    pdb_smiles_dict={}
    cleanedList = cleanlist(pdb_list)
    if len(cleanedList)!=0:
        ligs_chunks = [cleanedList[i:i + n] for i in range(0, len(cleanedList), n)]
        pdb_smiles_valid = reduce(list.__add__, [pdb_to_smiles( chunk, mock) for chunk in tqdm(ligs_chunks)])
        for i in pdb_smiles_valid:
            pdb=i["pdb_ligand"]
            smiles=i["smiles"]
            pdb_smiles.append(smiles)
            pdb_smiles_dict[pdb]=smiles
    else:
        pdb_smiles=[]
        pdb_smiles_dict={}
    return pdb_smiles_dict, pdb_smiles #Diccionario pdb:smiles    
    

def chembl_to_smiles(chembl_id, chembl_db):
    engine = create_engine(chembl_db)
    CHEMBL_VERSION = 28
    find_smilesbychembl = (f'''SELECT c.chembl_id, e.canonical_smiles
    FROM MOLECULE_DICTIONARY c 
    JOIN COMPOUND_STRUCTURES e ON c.MOLREGNO = e.MOLREGNO
    WHERE c.chembl_id= "{chembl_id}" ''')
    df_mol = pd.read_sql(find_smilesbychembl, engine)
    return df_mol

def chembl_smiles_dic(chembl_list, chembl_db):
    appended_data = []
    
    chembl_smiles_dict={}
    chembl_list=cleanlist(chembl_list)
    if len(chembl_list)!=0:
        for i in chembl_list:
            appended_data.append(chembl_to_smiles(i, chembl_db))
        df_chembl_to_smiles = pd.concat(appended_data, ignore_index=True)
        chembl_smiles=list(df_chembl_to_smiles["canonical_smiles"])
        chembl_smiles_dict=df_chembl_to_smiles.set_index('chembl_id').T.to_dict('records')
        for s in chembl_smiles_dict:
            chembl_smiles_dict=s
    else:
        chembl_smiles = []
        chembl_smiles_dict={}
    return chembl_smiles_dict, chembl_smiles

def smiles_filter(smiles_final, heavy_atoms=7): 
    '''Se eliminan todos los que presentan menos de 7 heavy atoms (C N O etc, no H)'''
    filtrado=[]
    rene=[]
    rene_dict={}
    filtrado_dict={}
    for mol, smiles in smiles_final.items():
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if not m:
            continue
        else:
            if m.GetNumHeavyAtoms()>=7:
                Chem.MolToSmiles(m)
                filtrado.append((mol, Chem.MolToSmiles(m)))
                filtrado_dict=dict((x, y) for x, y in filtrado)
            else:
                rene.append((mol, Chem.MolToSmiles(m)))
                rene_dict=dict((x, y) for x, y in rene)
    return filtrado_dict, rene_dict

def process_pdb_list(pdb_list_path):    
    #dn_chembl_assay_trusted
    pdb_list=[]
    pdb_list_file=open(pdb_list_path, "r")
    for pdb in pdb_list_file:
        pdb_list.append(pdb.strip())
    return pdb_list

def process_chembl_list(chembl_list_path):    
    #dn_chembl_assay_trusted
    chembl_list=[]
    chembl_list_file=open(chembl_list_path, "r")
    for chembl in chembl_list_file:
        chembl_list.append(chembl.strip())
    return chembl_list

def smiles_list_maker(path, chembl_db, n, mock=False):
    dire=open_directory(path)
    for i in tqdm(dire):
        sys.stderr.write(f'Running Uniprot ID= {i}\n')
        #print(i)
        filepath = path + os.sep + i
        if os.path.isfile(filepath+ os.sep +'smiles_dict_removed'):
            print(f'{i} already have smiles_list')
            pass
        else:
            pdb_list=process_pdb_list(filepath+ os.sep + 'pdb_lista_final.lst')
            pdb_smiles_dict, pdb_smiles=pdb_smiles_dic(pdb_list, n, mock)  
            chembl_list=process_chembl_list(filepath+ os.sep + 'chembl_lista_final.lst')
            chembl_smiles_dict, chembl_smiles=chembl_smiles_dic(chembl_list, chembl_db)
            smiles_dict_final = {}
            smiles_dict_final.update(pdb_smiles_dict)
            smiles_dict_final.update(chembl_smiles_dict)
            smiles_dict_filtered, smiles_dict_removed=smiles_filter(smiles_dict_final, 7)
            smiles_final=set(chembl_smiles+pdb_smiles)
            with open(filepath+ os.sep +'smiles_dict_final',"w") as h: 
                h.write(str(smiles_dict_final))
            with open(filepath+ os.sep +'smiles_list_final.lst',"w") as h: 
                h.write('\n'.join(smiles_final))
            with open(filepath+ os.sep +'smiles_dict_filtered',"w") as h: 
                h.write(str(smiles_dict_filtered))  
            with open(filepath+ os.sep +'smiles_dict_removed',"w") as h: 
                h.write(str(smiles_dict_removed))  


def main(mock=False):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help= 'Directory path', required=True)
    parser.add_argument('-db', '--dataset', help='Path to the directory of ChEMBL DB', default='chembl_28.db')
    parser.add_argument('-o_dict', '--output_dict', help= 'To smiles_dict_final name')
    parser.add_argument('-o_list', '--output_list', help= 'To change smiles_list_final.lst name')
    parser.add_argument('-o_dfiltered', '--output_dfiltered', help= 'To smiles_dict_filtered name')
    parser.add_argument('-o_dremoved', '--output_dremoved', help= 'To smiles_dict_removed name')
    parser.add_argument('-c_size', '--chunk_size', help= 'Set pdb chunk size.', default=50)
    parser.add_argument('-m', '--mock', help= 'To test without API query use mock=True', action='store_true')
    args=parser.parse_args()
    
    #chembl_db = 'sqlite:///'+ args.chembl_db
    if os.path.exists(args.dataset):
        chembl_db = 'sqlite:///'+os.path.abspath(args.dataset)
        sys.stderr.write(f'Running with chembl_db={chembl_db} and mock={args.mock}\n')
        if args.mock:
            try:
                test=smiles_list_maker(args.input, chembl_db, args.chunk_size, args.mock)
            except Exception as ex:
                #sys.stderr.write(str(ex))   
                print('Error processing ligands')
        else:
            list_maker=smiles_list_maker(args.input, chembl_db, args.chunk_size, args.mock)
    else:
        sys.stderr.write(f'{args.dataset} does not exist')


if __name__=='__main__':
    main()
