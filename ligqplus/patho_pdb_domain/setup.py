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
