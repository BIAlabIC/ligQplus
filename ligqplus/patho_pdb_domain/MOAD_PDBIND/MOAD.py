import json 
import argparse
import toMolar

def moad_parse(input_file):
    """Toma la base de datos MOAD y la parsea en un JSON con jerarquias Ligand->PDB->Residues."""
    input_file=open(input_file,"r")
    input_file_lines=input_file.readlines()
    input_file.close()

    compound_dict={}

    for line in input_file_lines:
        line_split=line.split(",")

        if all(value == "" or value == "\n" for value in line_split):
            continue

        if line_split[0]!="":
            continue

        if line_split[2]!="":
            pdb=line_split[2]
            continue

        compound=line_split[3].split(":")[0]

        chain=line_split[3].split(":")[1]
        resid=line_split[3].split(":")[2]
        status=line_split[4]
        if line_split[7] !="":
            afinity_standard=line_split[7]+line_split[8]
            afinity = toMolar.toMolar(float(line_split[7]),line_split[8])
            type_afinity=line_split[5]
            standard_relation=line_split[6]
        else:
            afinity="None"
            afinity_standard="None"
            type_afinity="None"
            standard_relation="None"

        if compound not in compound_dict.keys():
            compound_dict[compound]= { "pdbs": [] }
            record = { "name": pdb, "residues": []  }
            if type_afinity=="None":
                residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, "standard_relation": standard_relation}
            else:
                residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, type_afinity:afinity_standard, "standard_relation": standard_relation}
            record["residues"].append(residues)
            compound_dict[compound]["pdbs"].append(record)

        else:
            check=0
            for element in compound_dict[compound]["pdbs"]:
                if pdb == element["name"]:
                    if type_afinity=="None":
                        residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, "standard_relation": standard_relation}
                    else:
                        residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, type_afinity:afinity_standard, "standard_relation": standard_relation}
                    element["residues"].append(residues)
                    check=1
            if check==0:
                if type_afinity=="None":
                    residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, "standard_relation": standard_relation}
                else:
                    residues = {"chain": chain, "resid" : resid, "status": status, "standard_value": afinity, type_afinity:afinity_standard, "standard_relation": standard_relation}
                record = { "name": pdb, "residues": []  }
                record["residues"].append(residues)
                compound_dict[compound]["pdbs"].append(record)

    print(json.dumps(compound_dict, indent=4, sort_keys=True))
    return 0

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create Ligand->PDB->Residues structure json database from MOAD file')
    parser.add_argument("-i", '--moad_file', default='every_bind.csv', help="Provide input moad file and get json file")
    return parser

def main():
    parser=parse_arguments()
    args=parser.parse_args()
    moad_parse(args.moad_file)
    return 0

if __name__=='__main__':
    main()



