# ligQplus
Massive compound screening for pathogenic proteomes.

<img src="https://docs.google.com/drawings/d/1PApRMKCJE-YwFVnGwfm4VbXhjwNuiiGMwmvLRIKXxkU/export/png">

#MAPPING PDBs WITH PFAM IDs
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt
#LIST OF PDB ENTRIES WITH CO-CRISTALYZATE LIGANDS
http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig_pairs.lst
# DOWNLOAD MOAD DATABASE
http://bindingmoad.org/files/csv/every_bind.csv


# Web Service Client

#To help users get started with using the updated ChEMBL web services the existing web service client has also #been released. This is written in the Python programming language and is available to install from Python #Package Index by typing:
#chembl_version=27

chembl_webresource_client

#The client code is open and hosted on GitHub: https://github.com/chembl/chembl_webresource_client.

Chembl DB 

ChEMBL Database downloads, which includes Oracle, MySQL and PostgreSQL versions of the ChEMBL database, as well as SDF, FASTA and release note files. You will need it to search by pfam ID

ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/



mol_trg.py: returns all target related to the search compound according to API. INPUT=.txt list of Chemble_id (compound). OUTPUT= .csv table containing 3 columns trg_chembl_ID, 'pchembl_median', 'activity_comment'.

trg_mol.py: returns all compounds related to the search target according to API. INPUT= .txt list of Chemble_id (target). OUTPUT= .csv  table containing 4 columns mol_chembl_ID, 'pchembl_median', 'activity_comment', 'target_organism'.

tanimoto.py: returns all compounds related to the search SMILES according to the Tanimoto similarity score given (60% Default). INPUT= .txt SMILES list. OUTPUT= .csv  table containing 3 columns 'mol_chembl_id', 'similarity','smiles'

pfam_mol_assay.py: returns all compounds related to the search pfam ID according to ChEMBL data base via assay related information. INPUT= .txt PFAM ID list. OUTPUT= .csv  table containing 2 columns of 'mol_chembl_id': "Dudoso", "Confiable".

pfam_mol_mech.py: returns all compounds related to the search pfam ID according to ChEMBL data base via mechanism related information. INPUT= .txt PFAM ID list. OUTPUT= .csv  table containing 2 columns of 'mol_chembl_id': "Dudoso", "Confiable".

SMILES_file: SMILES list

OUT: File with two columns. First column ChEMBL_ID, second column SMILES. Separated by "\t". xPfamID + xTanimoto
<img src="https://docs.google.com/drawings/d/e/2PACX-1vSSwg9kpBGrZ5d2lJAgvReRPHrV0O1JAkZ2C8Mu9ui4F2FxBriT6iRT8mE1QZaTFPWPx9qbpNCMPNRf/pub?w=960&amp;h=720">
