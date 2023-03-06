from distutils.core import setup
setup(name='patho_pdb_domain',
		version='0.0.23',
		py_modules=['MOAD_PDBIND','extracts','ligand_from_pfam',
			    'ligand_from_pfam.domain_pdb_ligand','ligand_from_pfam.request_ligand_from_PDBe',
			    'extracts.extract_ligand_from_pdb','extracts.extract_pdb_from_domain',
			    'extracts.protein_id_extract_to_uniprot', 'MOAD_PDBIND.filter_MOAD'],
		scripts=['MOAD_PDBIND/MOAD.py','MOAD_PDBIND/filter_MOAD.py', 'MOAD_PDBIND/toMolar.py', 'MOAD_PDBIND/PDBBIND.py', 'extracts/extract_ligand_from_pdb.py'
			 ,'extracts/extract_pdb_from_domain.py','extracts/protein_id_extract_to_uniprot.py'
			 ,'ligand_from_pfam/domain_pdb_ligand.py','ligand_from_pfam/request_ligand_from_PDBe.py'],

		requires=['requests','argparse','biopython'],

		author='Federico Serral',
		license='MIT license',
		author_email='fedeserral92@gmail.com',
		description='Extract true ligands',
		url='https://github.com/fedeserral/patho_pdb_domain',
		long_description='',
		)



#from distutils.core import setup
import setuptools

setuptools.setup(name='target_chembl',
		version='0.0.7',
		scripts=['patho_chembl/chembldb_pfam_mech.py', 'patho_chembl/chembldb_pfam_assay.py', 'patho_chembl/mol_trg.py', 'patho_chembl/pfam_df_update.py', 'patho_chembl/pfam_mol_assay.py',
    		'patho_chembl/pfam_mol_mech.py', 'patho_chembl/pfam_trg_sql_assay.py', 'patho_chembl/pfam_trg_sql_mecanism.py', 'patho_chembl/tanimoto.py', 'patho_chembl/trg_mol.py', 'patho_chembl/trg_mol_funcion.py'],

		requires=['requests','argparse', 'chembl_webresource_client', 'pandas'],

		author='Florencia A. Castello',
		license='MIT license',
		author_email='florencia.castelloz@gmail.com',
		description='Simple interface for ChEMBL DB',
		url='https://github.com/florenciacastello/target_chembl',
      		packages=setuptools.find_packages(),
		long_description='',
		python_requires='>=3.6'
		)
