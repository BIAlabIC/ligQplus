# ligQplus
Massive compound screening for pathogenic proteomes.

## Work scheme

  In order to carry out this pipeline, the first thing to do is create a conda environment where you can install a series of necessary packages. See here more information about [Anaconda](https://www.anaconda.com/download/). In addition, it is necessary to have a version of python 3 at least.
  
### Download and pre-requisites 
```
$ git clone https://github.com/BIAlabIC/ligQplus.git
```
### Create execution environment:
```
$ conda create --name ligQ_plus
```
### Activate environment
```
$ conda activate ligQ_plus
```
Once in the ligQ plus environment install:
```
$ pip install tqdm
$ pip install pandas
$ pip install Bio
$ pip install chembl_webresource_client

```
In order for your Linux system to find the given command/executable file, we need to make sure to run the following command in the using terminal before starting work:
```
$ export PYTHONPATH=$PYTHONPATH:/path/to/ligqplus
$ export PATH=$PATH:/path/to/ligqplus
```

## Create databases
 In addition to the scripts, it is necessary to have a series of Databases. Some scripts search using databases locally via SQL such as ChEMBL and others via URL/API. Some databases should be previously generated to save computing time. It is necessary to have the space to store these databases.
### MOAD
The every_bind.csv file is generated based on experimental assays and contains curated information about interactions between ligands and proteins available in the PDB database.
```
$ wget -O moad.csv "http://bindingmoad.org/files/csv/every_bind.csv"

$ python3 /patho_pdb_domain/MOAD_PDBIND/MOAD.py > moad.json
```
### CHEMBL
The [ChEMBL](https://www.ebi.ac.uk/chembl/) database can be at least 20gb in size, please check if you have enough space before downloading.
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/#last_version#/
```
From this base, you must create 2 tables that summarize the information of [PFAM](http://pfam.xfam.org/) domains per target molecule and the ligands that each presents:
```
$ python3 target_chembl/patho_chembl/pfam_trg_sql_assay.py -db chembl_#.db > pfam_assay_##.csv

$ python3 target_chembl/patho_chembl/pfam_trg_sql_mech.py -db chembl_#.db > pfam_mech_#.csv
```
### Mapping pfam-pdb and pdb-ligands
The pdb_pfam_mapping.txt file contains information about the PFAM domains (start/end positions) of each of the protein structures deposited in the PDB.

The lig_pairs.lst file contains information about the ligands co-crystallized with each of the protein structures deposited in the PDB.
```
$ wget -O pdb_pfam_mapping.txt "http://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt"

$ wget -O lig_pairs.lst "http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig_pairs.lst"
```

### Prepare proteome and extract Uniprot IDs
> Note: "Proteome_uniprot.fasta" refers to all proteins in uniprot DataBase related to the organism of interest. Must be download from [Uniprot](https://www.uniprot.org)
```
$ makeblastdb -in Proteome_uniprot.fasta -dbtype prot

$ blastp -db Proteome_uniprot.fasta -query organism.fasta -qcov_hsp_perc 80 -num_threads 4 -max_hsps 1 -outfmt 6 > file.txt

$ tr '.' ',' < file.txt | awk '$3>=90' > organism_filter.txt (filter by identity)

$ awk -F " " '{print $1, "", $2}' organism_filter.txt > organism_final.txt (keep the first 2 columns)
```
> 'organism_final.txt' It will finally have 2 columns. The first containing the locus tag of each protein of the organism of interest and the second column with the respective UniprotID. There can be multiple Uniprot IDs for each locus tag.

## Scripts:

> note: Make sure to have all the files, DataBases and scripts in the same folder before proceeding. Including: organism/XXX_final.txt, pfam_assay_##.csv, pfam_mech_##.csv, chembl_##.db, moad.json, etc.

### Run

1. prepocessing_proteome_liqQ.py
```
  $ python3 prepocessing_proteome_liqQ.py -i /path/to/organism/organism_final.txt -o organism -n organism
```
2. Create organism directory (mkdir organism inside path/to/organism/)

  * make_dic.py (create directories with the list of Uniprot IDs to test)
```
  $ python3 make_dic.py -i organism_lt.txt (output from 1.) -o organism/organism/
```
3. proteome_compound.py (returns a list of ligands by Uniprot ID)
```
  $ python3 proteome_compound.py -i organism_lt.txt (output from 1.) -o organism/organism/
```
4. analysis_result.py (Takes the previous output and combines the results by ChEMBL and PDB)
```
  $ python3 analysis_result.py -i organism/organism/
```
5. smile_list.py (takes the previous output and converts the ligands into their SMILES. It requires chembl_#.db for SQL)
```
  $ python3 smile_list.py -i organism/organism/ 
```
6. lt_csv_maker/lt_id_smiles/lt_smiles_maker.py (creates the final outputs by grouping the ligands by lt)
```
  $ python3 lt_csv_maker.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_ligands.csv

  $ python3 lt_id_smiles.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_id_smiles.csv

  $ python3 lt_smiles_maker.py -i organism/organism/ -dict organism/organism_final.json (output from prepocessing_proteome_liqQ.py) -o organism_lt_smiles.csv
```
7. clustering.py (Clusters the ligands by lt and searches for commercially available ligands in the Zinc database.)
```
  $ python3 clustering.py -i organismo_lt_id_smiles.csv -o organism/organism/ 

```
