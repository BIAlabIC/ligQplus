#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 15:17:22 2021

@author: fleer
"""
import pandas as pd
import requests
from bioservices import *
from bioservices import UniChem
import os
import argparse
from tqdm import tqdm
import sys
from functools import reduce
from requests.adapters import HTTPAdapter 
from requests.packages.urllib3.util.retry import Retry
from requests.structures import CaseInsensitiveDict


#from io import StringIO

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
    
def process_pdb_ligand(pdb_ligand_path):
    #PDB_ligand
    pdb_ligand_final=[]
    with open(pdb_ligand_path) as pdb:
        df_pdb_ligands=pd.read_table(pdb)
        if "ligand" in df_pdb_ligands.columns:
            pdb_ligand=list(df_pdb_ligands["ligand"])
            pdb_ligand=cleanlist(pdb_ligand)
            for ligand in pdb_ligand:
                ligand=str(ligand)
                if '-' not in ligand:
                    pdb_ligand_final.append(ligand)
            pdb_ligand_final=set(pdb_ligand_final)
        else:
            pdb_ligand_final=set([])
    return pdb_ligand_final

def process_pdb_ligand_valid(pdb_ligand_valid_path):
    #PDB_ligand_valid
    pdb_ligand_valid_final=[]
    if os.path.getsize(pdb_ligand_valid_path) != 0:
        with open(pdb_ligand_valid_path) as pdb_valid:
            df_pdb_ligand_valid=pd.read_table(pdb_valid, sep='\s+', header=None)
            if df_pdb_ligand_valid.empty:
                pdb_ligand_valid_final=set([])
            else:
                pdb_ligand_valid=list(df_pdb_ligand_valid[0])
                pdb_ligand_valid=cleanlist(pdb_ligand_valid)
                for ligand in pdb_ligand_valid:
                    ligand=str(ligand)
                    if '-' not in ligand and len(ligand)!=1:
                        pdb_ligand_valid_final.append(ligand)
                pdb_ligand_valid_final=set(pdb_ligand_valid_final)
    else:
        pdb_ligand_valid_final=set([])
    return pdb_ligand_valid_final

def process_dn_pdb_ligand(dn_pdb_ligand_path):
    #dn_PDB_ligand
    dn_pdb_ligand_final=[]
    with open(dn_pdb_ligand_path) as dn_pdb:
        df_dn_pdb_ligands=pd.read_table(dn_pdb)
        df_dn_pdb_ligands=df_dn_pdb_ligands.dropna()
        if "ligand" in df_dn_pdb_ligands.columns:
            dn_pdb_ligand=list(df_dn_pdb_ligands["ligand"])
            dn_pdb_ligand=cleanlist(dn_pdb_ligand)
            for ligand in dn_pdb_ligand:
                if '-' not in ligand:
                    dn_pdb_ligand_final.append(ligand)
            dn_pdb_ligand_final=set(dn_pdb_ligand_final)
        else:
            dn_pdb_ligand_final=set([])
    return dn_pdb_ligand_final

def process_dn_pdb_ligand_valid(dn_pdb_ligand_valid_path):
    #dn_PDB_ligand_valid
    dn_pdb_ligand_valid=[]
    dn_pdb_ligand_valid_final=[]
    if os.path.getsize(dn_pdb_ligand_valid_path) != 0:
        with open(dn_pdb_ligand_valid_path) as dn_pdb_valid:
            df_dn_pdb_ligands_val=pd.read_table(dn_pdb_valid, sep='\s+', names=["Ligand", "PDB"])
            df_dn_pdb_ligands_val=df_dn_pdb_ligands_val.dropna()
            if df_dn_pdb_ligands_val.empty:
                dn_pdb_ligand_valid=set([])
            else:
                dn_pdb_ligand_valid=list(df_dn_pdb_ligands_val["Ligand"])
                dn_pdb_ligand_valid=cleanlist(dn_pdb_ligand_valid)
                for ligand in dn_pdb_ligand_valid:
                    ligand=str(ligand)
                    if '-' not in ligand and len(ligand)!=1:
                        dn_pdb_ligand_valid_final.append(ligand)
            dn_pdb_ligand_valid_final=set(dn_pdb_ligand_valid_final)
    return dn_pdb_ligand_valid_final

def process_chembl_comp(chembl_comp_path):
    #chembl_comp
    chembl_comp=[]
    with open(chembl_comp_path) as chembl:
        df_chembl_comp=pd.read_table(chembl)
        if "molecule_chembl_id" in df_chembl_comp.columns:
            chembl_comp=set(list(df_chembl_comp["molecule_chembl_id"]))
        else:
            chembl_comp=[]
    return chembl_comp

def process_dn_chembl_assay_t(dn_chembl_assay_t_path):    
    #dn_chembl_assay_trusted
    dn_chembl_assay_t=[]
    if os.path.getsize(dn_chembl_assay_t_path) != 0:
        dn_chembl_assay_t=[]
        dn_chembl_a=open(dn_chembl_assay_t_path, "r")
        for chembl in dn_chembl_a:
            dn_chembl_assay_t.append(chembl.strip())
    else:
        dn_chembl_assay_t=[]
    return dn_chembl_assay_t

def process_dn_chembl_assay_u(dn_chembl_assay_u_path):
    #dn_chembl_assay_uncertain
    dn_chembl_assay_u=[]
    if os.path.getsize(dn_chembl_assay_u_path) != 0:
        dn_chembl_a=open(dn_chembl_assay_u_path, "r")
        for chembl in dn_chembl_a:
            dn_chembl_assay_u.append(chembl.strip())
    else:
        dn_chembl_assay_u=[]
    return dn_chembl_assay_u

def process_dn_chembl_mech_t(dn_chembl_mech_t_path):
    #dn_chembl_mech_trusted
    dn_chembl_mec_t=[]
    if os.path.getsize(dn_chembl_mech_t_path) != 0:
        dn_chembl_m_t=open(dn_chembl_mech_t_path, 'r')
        for chembl in dn_chembl_m_t:
            dn_chembl_mec_t.append(chembl.strip())
    else:
         dn_chembl_mec_t=[]   
    return dn_chembl_mec_t

def process_dn_chembl_mech_u(dn_chembl_mech_u_path):
    #dn_chembl_mech_uncertain
    dn_chembl_mec_u=[]
    if os.path.getsize(dn_chembl_mech_u_path) != 0:
        dn_chembl_m_u=open(dn_chembl_mech_u_path, 'r')
        for chembl in dn_chembl_m_u:
            dn_chembl_mec_u.append(chembl.strip())
    else:
        dn_chembl_mec_u=[]
    return dn_chembl_mec_u

def pdb_ligand_data(ligands, mock=False): #usar el de batch que se pueden usar listas
    new_data = []
    if mock:
        return new_data
    retry_strategy = Retry(total=5, status_forcelist=[104, 429, 500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session() 
    http.mount("https://", HTTPAdapter(max_retries=retry_strategy))
    r = http.post("https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/",data=",".join(ligands))
    if r.ok:
        data =  r.json()  
        for k,v in data.items():
            r = v[0]            
            r2 = {"pdb_ligand":k}
            #if len(r["smiles"]):
            #    r2["smiles"] = r["smiles"][0]["name"]
            if "chembl_id" in r:
                r2["chembl_id"] = r["chembl_id"]            
                new_data.append(r2)
        return new_data
    else:
        #print(r)
        return new_data
    #raise Exception(r.text)

def pdb_to_chembl(pdb_list, n=50, mock=False):
    pdb_chembl=[]
    cleanedList = cleanlist(pdb_list)
    if len(cleanedList)!=0:
        ligs_chunks = [cleanedList[i:i + n] for i in range(0, len(cleanedList), n)]
        pdb_chembl_valid = reduce(list.__add__, [pdb_ligand_data( chunk ) for chunk in tqdm(ligs_chunks)])
        for data in pdb_chembl_valid:
            if data['chembl_id']!=None:
                pdb_chembl.append(data['chembl_id'])
    else:
        pdb_chembl=[]
    return pdb_chembl

def chembl_to_pdb(chembl_list):  
#    print(chembl_list)
    chembl_to_pdb=[]
    if len(chembl_list)!=0:
        for chembl_id in tqdm(chembl_list):
            data_chembl={"compound": chembl_id, "sourceID": 1,"type": "sourceID"}
            retry_strategy = Retry(total=3, status_forcelist=[429, 500, 502, 503, 504], 
                       method_whitelist=["HEAD", "GET", "OPTIONS", "POST"])
            adapter = HTTPAdapter(max_retries=retry_strategy)
            http = requests.Session() 
            http.mount("https://", adapter) 
            http.mount("http://", adapter)
            headers = {
                'accept': 'application/json',
                # Already added when you pass json= but not when you pass data=
                # 'Content-Type': 'application/json',
            }
            r = http.post("https://www.ebi.ac.uk/unichem/api/v1/compounds", headers=headers, json=data_chembl)
            if r.ok:
                data =  r.json()  
                for k,v in data.items():
                    if k == "compounds":
                        #print(v)
                        if v:
                            compound_inf=v[0]
                            for info, info_list in compound_inf.items():
                                if info=='sources':
                                    # print(info_list)
                                    for source_info in info_list:
                                        # for source_k, source_v in source_info.items():
                                        if source_info['id']==3:
                                            pdb_ligandid=source_info['compoundId']
                                            chembl_to_pdb.append(pdb_ligandid)
            # print(chembl_to_pdb)         
            # raise Exception(r.text)

    return chembl_to_pdb
    


def dif_pdb_chembl(pdb_valid_sum, chembl_to_pdb): #Difference between pdb_ID list and pdb_id from chembl_id.
    dif_pdb_chembl_final=[]
    dif_pdb=[]
    pdbs=set(pdb_valid_sum)
    chembl_pdb=set(chembl_to_pdb)
    diference=pdbs.difference(chembl_pdb)
    dif_pdb.append(diference)
    if len(dif_pdb):
        dif_pdb_chembl_final=[pdb for x in dif_pdb for pdb in x]
        dif_pdb_chembl_final=cleanlist(dif_pdb_chembl_final)
    else:
        dif_pdb_chembl_final=[]
    return dif_pdb_chembl_final

def dif_chembl_pdb(chembl, dif_pdb_chembl_final_to_chembl):
    difs_chembl=[]
    difs_chembl_final=[]
    pdbs_to_chembl=set(dif_pdb_chembl_final_to_chembl)
    difs_chembl.append(chembl.difference(pdbs_to_chembl))
    if len(difs_chembl)!=0:
        difs_chembl_final=[chemb for x in difs_chembl for chemb in x]
        difs_chembl_final=cleanlist(difs_chembl_final)
    else:
        difs_chembl_final=[]
    return difs_chembl_final


def final_list_maker(path, n, mock=False):
    dire=open_directory(path)
    for i in tqdm(dire):
        sys.stderr.write(f'Running Uniprot ID= {i}\n')
        filepath = path + os.sep + i
        if os.path.isfile(filepath+ os.sep +'chembl_lista_final.lst'):
            print(f'{i} already have final_list')
            pass
        else:
            pdb_ligands=process_pdb_ligand(filepath+ os.sep + 'pdb_ligands.tbl')
            pdb_ligand_valid=process_pdb_ligand_valid(filepath+ os.sep + 'pdb_ligands_valid.tbl')
            dn_pdb_ligand=process_dn_pdb_ligand(filepath+ os.sep + 'dn_ligands.tbl')
            dn_pdb_ligand_valid=process_dn_pdb_ligand_valid(filepath+ os.sep + 'dn_ligands_valid.lst')
            chembl_comp=process_chembl_comp(filepath+ os.sep + 'chembl_comp.tbl')
            dn_chembl_assay_t=process_dn_chembl_assay_t(filepath+ os.sep + 'dn_chembl_assay_trusted.lst')
            dn_chembl_assay_u=process_dn_chembl_assay_u(filepath+ os.sep + 'dn_chembl_assay_uncertain.lst')
            dn_chembl_mec_t=process_dn_chembl_mech_t(filepath+ os.sep + 'dn_chembl_mec_trusted.lst')
            dn_chembl_mec_u=process_dn_chembl_mech_u(filepath+ os.sep + 'dn_chembl_mec_uncertain.lst')
            
            #lists
            chembl=set(dn_chembl_mec_t+dn_chembl_assay_t+list(chembl_comp))
            chembl_list=list(chembl)
            chembl_list=cleanlist(chembl_list)
            pdb=list(dn_pdb_ligand_valid)+list(pdb_ligand_valid)
            pdb=cleanlist(pdb)
            
            #From Chembl_ID to PDB_ID
            chembl_to_pdb_list=chembl_to_pdb(chembl_list)
            
            #PDB Final List
            dif_pdb_chembl_final=dif_pdb_chembl(pdb, chembl_to_pdb_list)
            
            #PDb final to Chembl
            dif_pdb_chembl_final_to_chembl=pdb_to_chembl(dif_pdb_chembl_final, n, mock)
                        
            #Chembl final list
            dif_chembl_pdb_final=dif_chembl_pdb(chembl, dif_pdb_chembl_final_to_chembl)
            
            #Final list
            final_list=dif_pdb_chembl_final+dif_chembl_pdb_final
            final_list = cleanlist(final_list)
            with open(filepath+ os.sep +'final_list.lst',"w") as h: 
                h.write('\n'.join(final_list))
            with open(filepath+ os.sep +'pdb_lista_final.lst',"w") as h: 
                h.write('\n'.join(dif_pdb_chembl_final))        
            with open(filepath+ os.sep +'chembl_lista_final.lst',"w") as h: 
                h.write('\n'.join(dif_chembl_pdb_final)) 
            

def main(mock=False):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help= 'Directory path', required=True)
    parser.add_argument('-o_pdb', '--output_pdb', help= 'To change pdb_lista_final.lst name')
    parser.add_argument('-o_chembl', '--output_chembl', help= 'To change chembl_lista_final.lst name')
    parser.add_argument('-o_final', '--output_final', help= 'To change final_list.lst name')
    parser.add_argument('-c_size', '--chunk_size', help= 'Set pdb chunk size.', default=50)
    parser.add_argument('-m', '--mock', help= 'To test without API query use mock=True', action='store_true')
    args=parser.parse_args()
    if args.mock:
        test=final_list_maker(args.input, args.chunk_size, args.mock)
    else:
        list_maker=final_list_maker(args.input, args.chunk_size, args.mock)


if __name__=='__main__':
    main()
