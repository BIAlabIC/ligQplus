import argparse
from collections import defaultdict
import json
import toMolar

def pdb_bind(input_file):
    mapping= open (input_file,"r")
    mapping_lines=mapping.readlines()
    mapping.close()
    db_pdb_bind= defaultdict(list)
    for line in mapping_lines:
        PDB , resolution, year, log_assay, assay, reference, raw_ligand = line.strip().split(" ")[:7]
        standard_relation= [x for x in "=<>~" if x in assay][0]
        medida, valor = assay.split(standard_relation)
        ligand= raw_ligand.replace("_", "")

        if "-mer" not in ligand:
            multiple_ligands= [x for x in ["-","/","&"] if x in ligand]
            if multiple_ligands:
                multiple_ligands = multiple_ligands[0]
                ligands= ligand.split(multiple_ligands)

            else:
                multiple_ligands = ""
                ligands= [ligand]

            for compound in ligands:

                if compound not in db_pdb_bind:
                    db_pdb_bind[compound] = {  "pdbs": [ ] }

                #Convierto las unidades a Molar
                valor_valor=float(valor[:-2])
                valor_unidad=valor[-2:]
                valor_Molar=toMolar.toMolar(valor_valor,valor_unidad)

                record =  {  "name": PDB, "affinity": { medida:valor, "standard_relation":standard_relation, "standard_value":valor_Molar }  }
                if multiple_ligands:
                    record["co_ligand"]={"relationship":multiple_ligands, "ligand":
                    (ligands[0] if ligands[0] != compound else ligands[1])}
                db_pdb_bind[compound] ["pdbs"].append( record )

    print(json.dumps(db_pdb_bind, indent=4, sort_keys=True))
    return 0

pdb_bind("pdb_bind_db.2019")



