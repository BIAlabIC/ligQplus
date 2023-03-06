#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 09:36:39 2021

@author: fleer
"""
from tqdm.auto import tqdm
import pandas as pd
import uniprot as uni
from patho_chembl.trg_mol_funcion import trg_to_mol
from ligand_from_pfam.domain_pdb_ligand import ligands_from_domain
from extracts.extract_ligand_from_pdb import ligands_from_pdb, pdb_ligands_mapping
from MOAD_PDBIND.filter_MOAD import filter_ligands
from ligand_from_pfam.request_ligand_from_PDBe import pdb_ligand_data_batch
from patho_chembl.pfam_mol_assay import fix_dataset, search_bypfam
from patho_chembl.pfam_mol_mech import search_bypfam as search_bypfam_mech
import json
import argparse
import os
from requests.exceptions import RequestException 



def proteome_compound(uniprot_list, path_dic, db_assay, db_mech, pdb_pfam_dataset, lig_pairs, moad):
    db_dn_assay = pd.read_csv(db_assay) 
    db_dn_mech = pd.read_csv(db_mech)
    fix_db_assay = fix_dataset(db_dn_assay)
    fix_db_mech = fix_dataset(db_dn_mech)
    #pdb_pfam=open(pdb_pfam_dataset)
    #lig_pairs=open(lig_pairs)
    with open (lig_pairs) as lig_pairs:
        pdb_ligand_mapping_dict = pdb_ligands_mapping(lig_pairs)
    with open(moad) as j: 
        moad_json = json.load(j)
    for input_unip in tqdm(uniprot_list):
        input_unip=input_unip.strip()
        path = path_dic+ input_unip
        path_empty=os.listdir(path)
        print('Running ', input_unip)
        if len(path_empty) == 0:
            pdbs = uni.map(input_unip,f="UniProtKB_AC-ID",t='PDB')
            print(input_unip, pdbs)
            pdb_list=[]
            if "failedIds" not in pdbs.keys():
                for k in pdbs.keys():
                    for v in pdbs.values():
                        for element in v:
                            print(element)
                            res = element.get('to')
                            pdb_list.append(res)
            #print (res)

            chembl = uni.map(input_unip,"UniProtKB_AC-ID","ChEMBL")
            #print(chembl)                
            chembl_comp =[] 
            chembl_list=[]
            if "failedIds" not in chembl.keys():
                for k in chembl.keys():
                    for v in chembl.values():
                        for element in v:
                            chembl_list = element.get('to')
                            #print (res)
                            
            domains = uni.pfam_from_uniprot(input_unip)
            #print(domains)   
            if chembl_list:
                for chembl_target_id in chembl_list:
                    chembl_comp.append(trg_to_mol(chembl_target_id))
                chembl_comp=pd.concat(chembl_comp, ignore_index=True)
            else: 
                chembl_comp=pd.DataFrame()
                
            if pdb_list:
                #print(pdbs)
                pdb_ligands = ligands_from_pdb(pdb_list,pdb_ligand_mapping_dict)         
                pdb_ligands_valid = filter_ligands(pdb_ligands .to_records(index=False),moad_json)
            else:
                pdb_ligands =pd.DataFrame()
                pdb_ligands_valid =pd.DataFrame()
            with open (pdb_pfam_dataset) as pdb_pfam:
                dn_ligands = ligands_from_domain(domains, pdb_pfam)
            dn_ligands_valid = filter_ligands(dn_ligands[["ligand","pdb"]].to_records(index=False),moad_json)
    
            dn_chembl = search_bypfam(domains, fix_db_assay) #busca por PFAM ID y devuelve los chemblID de los compuestos relacionados con los target que poseen esos pfamID 
            dn_mech = search_bypfam_mech(domains, fix_db_mech)
    
            with open(f'{path}/chembl_comp.tbl',"w") as h:
                chembl_comp.to_csv(h,sep="\t")
            with open(f'{path}/pdb_ligands.tbl',"w") as h:
                pdb_ligands.to_csv(h,sep="\t")
            with open(f'{path}/pdb_ligands_valid.tbl',"w") as h:
                h.write('\n'.join('{}    {}'.format(ligand[0],ligand[1]) for ligand in pdb_ligands_valid))
            with open(f'{path}/dn_ligands.tbl',"w") as h:
                dn_ligands.to_csv(h,sep="\t")
            with open(f'{path}/dn_ligands_valid.lst',"w") as h: # imprimir en 2 columnas
                h.write('\n'.join('{}    {}'.format(ligand[0],ligand[1]) for ligand in dn_ligands_valid))    
           
            with open(f'{path}/dn_chembl_mec_trusted.lst',"w") as h: # separar los sets en 2 archivos
                h.write('\n'.join(x for x in dn_mech[0]))  #[x for x in dn_mech[0]]
            with open(f'{path}/dn_chembl_mec_uncertain.lst',"w") as h: # separar los sets en 2 archivos
                h.write('\n'.join(x for x in dn_mech[1])) 
            with open(f'{path}/dn_chembl_assay_trusted.lst',"w") as h: # separar los sets en 2 archivos
                h.write('\n'.join(x for x in dn_chembl[0])) 
            with open(f'{path}/dn_chembl_assay_uncertain.lst',"w") as h: # separar los sets en 2 archivos
                h.write('\n'.join(x for x in dn_chembl[1]))
        else:
            print('Already run')
            continue

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help='Input. File containing list of uniprot ID', type=argparse.FileType('r'), required=True, )
    parser.add_argument('-o', '--output', help= 'Directory path', required=True)
    parser.add_argument("-chembl_assay_db", "--chembl_assay_dataset", help="Path to the directory of ChEMBL DB", type=argparse.FileType('r'), default= "pfam_assay_28.csv")
    parser.add_argument("-chembl_mech_db", "--chembl_mech_dataset", help="Path to the directory of ChEMBL DB",  type=argparse.FileType('r'), default= "pfam_mech_28.csv")
    parser.add_argument("-pdb_pfam", "--pdb_pfam_mapping", help="Path to the directory of PDB Pfam mapping", default= "pdb_pfam_mapping.txt")
    parser.add_argument("-lig_pairs", "--lig_pairs_dataset", help="Path to the directory of lig pairs list", default= "lig_pairs.lst")
    parser.add_argument("-moad_db", "--moad_dataset", help="Path to the directory of moad DB", default= "moad.json")
        
    args=parser.parse_args()
    try:
        compounds = proteome_compound(args.input, args.output, args.chembl_assay_dataset, args.chembl_mech_dataset,
                                  args.pdb_pfam_mapping, args.lig_pairs_dataset, args.moad_dataset)
    except RequestException as e:
        print('Error using API. Try again.')
        raise
        


 
if __name__=='__main__':
    main()
