# ligQplus
Massive compound screening for pathogenic proteomes.

## Work scheme

### Create execution environment:
```
$ conda create --name ligQ_plus
```
### Activate environment
```
$ conda activate ligQ_plus
```
### Download and pre-requisites 
```
$ git clone https://github.com/BIAlabIC/ligQplus.git
$ export PYTHONPATH=$PYTHONPATH:/path/to/ligqplus
$ export PATH=$PATH:/path/to/ligqplus
```
## Create databases

### MOAD
```
$ wget -O moad.csv "http://bindingmoad.org/files/csv/every_bind.csv"

$ python /patho_pdb_domain/MOAD_PDBIND/MOAD.py > moad.json
```
### CHEMBL
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/#last_version#/

$ python target_chembl/patho_chembl/pfam_trg_sql_assay.py -db chembl_#.db > pfam_assay_##.csv

$ python target_chembl/patho_chembl/pfam_trg_sql_mech.py -db chembl_#.db > pfam_mech_#.csv
```
### Mapping pfam-pdb and pdb-ligands
```
$ wget -O pdb_pfam_mapping.txt "http://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt"

$ wget -O lig_pairs.lst "http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig_pairs.lst"
```

### Prepare proteome and extract Uniprot IDs
> Note: "Proteome_uniprot.fasta" refers to all proteins in uniprot DataBase related to the organism of interest. Must be download from [Uniprot](https://www.uniprot.org)
```
$ makeblastdb -in Proteome_uniprot.fasta -dbtype prot

$ blastp -db Proteome_uniprot.fasta -query XXX.fasta -qcov_hsp_perc 80 -num_threads 4 -max_hsps 1 -outfmt 6 > file.txt

$ tr '.' ',' < file.txt | awk '$3>=90' > XXX_filter.txt (filter by identity)

$ awk -F " " '{print $1, "", $2}' XXX_filter.txt > XXX_final.txt (keep the first 2 columns)
```

## Scripts:

> note: Make sure to have all the files and scripts in the same folder before proceeding. Including: organism/XXX_final.txt, pfam_assay_##.csv, pfam_mech_##.csv, chembl_##.db, moad.json, etc.

### Run

1. prepocessing_proteome_liqQ.py
```
  $ python prepocessing_proteome_liqQ.py -i /path/to/organism/XXX_final.txt -o organism -n organism
```
2. Create organism directory (mkdir organism inside path/to/organism/)

  * make_dic.py (create directories with the list of Uniprot IDs to test)
```
  $ python make_dic.py -i organism_lt.txt (output from 1.) -o organism/organism/
```
3. proteome_compound.py (returns a list of ligands by Uniprot ID)
```
  $ python proteome_compound.py -i organism_lt.txt (output from 1.) -o organism/organism/
```
4. analysis_result.py (Takes the previous output and combines the results by ChEMBL and PDB)
```
  $ python analysis_result.py -i organism/organism/
```
5. smile_list.py (takes the previous output and converts the ligands into their SMILES. It requires chembl_#.db for SQL)
```
  $ python smile_list.py -i organism/organism/ 
```
6. lt_csv_maker/lt_id_smiles/lt_smiles_maker.py (creates the final outputs by grouping the ligands by lt)
```
  $ python lt_csv_maker.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_ligands.csv

  $ python lt_id_smiles.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_id_smiles.csv

  $ python lt_smiles_maker.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_smiles.csv
```
7. clustering.py (Clusters the ligands by lt and searches for commercially available ligands in the Zinc database.)
```
  $ python clustering.py -i organismo_lt_id_smiles.csv -o organism/organism/ 

```
