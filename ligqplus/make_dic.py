# -*- coding: utf-8 -*-


import os, argparse, sys
def crear_directorio(path,lista_unip):
    for uniprot_id in lista_unip:
        uniprot_id= uniprot_id.strip()
        prot_dire =  path + uniprot_id
        isdir = os.path.isdir(prot_dire)
        if isdir:
            pass
        else:
            os.mkdir(prot_dire) 
        

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help='Input. File containing list of uniprot ID', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o', '--output', help= 'Directory path', required=True)
    args=parser.parse_args()
    prote_dire= crear_directorio(args.output, args.input)
 
if __name__=='__main__':
    main()
    
