# ligQplus

Massive compound screening for pathogenic proteomes.

# Work scheme

# Create execution environment:

conda create --name ligQ_plus

# Activate environment

conda activate ligQ_plus

# Download and pre-requisites 
git clone https://github.com/BIAlabIC/ligQplus.git
export PYTHONPATH=$PYTHONPATH:/path/to/ligqplus
export PATH=$PATH:/path/to/ligqplus

# Create databases

# MOAD
wget -O moad.csv "http://bindingmoad.org/files/csv/every_bind.csv"

python /patho_pdb_domain/MOAD_PDBIND/MOAD.py > moad.json
# CHEMBL
wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/#last_version#/

python target_chembl/patho_chembl/pfam_trg_sql_assay.py -db chembl_#.db > pfam_assay_##.csv

python target_chembl/patho_chembl/pfam_trg_sql_mech.py -db chembl_#.db > pfam_mech_#.csv

# Mapping pfam-pdb and pdb-ligands

wget -O pdb_pfam_mapping.txt "http://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt"

wget -O lig_pairs.lst "http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig_pairs.lst"


# Prepare proteome and extract Uniprot IDs

makeblastdb -in Proteome_uniprot.fasta -dbtype prot



























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
<img src="https://docs.google.com/drawings/d/e/2PACX-1vSSwg9kpBGrZ5d2lJAgvReRPHrV0O1JAkZ2C8Mu9ui4F2FxBriT6iRT8mE1QZaTFPWPx9qbpNCMPNRf/pub?

