#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse
import os
from collections import defaultdict

def lt_list_maker(blast_proteome_output, output_path, organism_name):
    protmap = defaultdict(list) 
    protmap2 = defaultdict(list)
    protmap3 = {}    
    path= blast_proteome_output
    for x in open(path): 
        prot,unip = x.split("  ")[:2]
        prot = prot.strip()
        unip = unip.strip()
        protmap[prot].append(unip) 
    
    for k,v in protmap.items(): 
        sp =[]
        for x in v:
            if x.startswith("sp"):
                sp.append(x)
                protmap2[k] = sp 
            else:
                vv =[]
                if len(v) == 1 :
                    vv.append(v[0])
                elif not sp:
                    vv.append(v[0])
                else:
                    vv=sp
            
            protmap2[k] = vv
             
    for key, value in protmap2.items():
        for v in value:
            v = v.split('|')[1]
        #print(v)
            protmap3[v] = key

##### dictionary {uniprot:[locus_tags],...}
    with open(f"{output_path}/{organism_name}_final.json","w") as h: 
        json.dump(protmap3,h)
    protmap3_keys= list(protmap3.keys())

#### list of locus tag split by \n
    for keys in protmap3_keys:
        with open(f"{output_path}/{organism_name}_lt.txt","a") as h: 
        	h.write(keys + "\n") 


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--blast_output", help=".txt file with information from blast between proteins and uniprot db",required = True)
    parser.add_argument("-o","--output_path", help="Directory path)",required=True)
    parser.add_argument("-n","--organism_name", help="name of specific organism proteome to run ",required = True)
    
    args = parser.parse_args()
    try:
        isExist_directory = os.path.exists(args.output_path)
        isExist_txt_file = os.path.exists(args.blast_output)
        if isExist_directory and isExist_txt_file:
            lt_list = lt_list_maker(args.blast_output, args.output_path, args.organism_name)
    except RequestException as e:
        print('Check path or input data and try again.')
        raise

    
if __name__=='__main__':
    main()
