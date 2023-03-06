import argparse
from chembl_webresource_client.new_client import new_client
from io import StringIO
import pandas as pd
from pandas import json_normalize

#ejemplo=['COC1=C2C(OC=C2)=NC3=CC(OC)=C(OC)C=C31', 'COC1=C2C(OC=C2)=NC3=C(OC)C4=C(OCO4)C=C31']

def chembl_similarity (handle, tanimoto_similarity):
	similarity=[]
	if handle:
		for smiles in handle:
			smiles = smiles.strip()
			similarity=similarity+list(new_client.similarity.filter(smiles = smiles, similarity = tanimoto_similarity).only('molecule_chembl_id',
			'similarity', 'molecule_structures'))
	return similarity

def Main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input', help='Input smiles_file', type=argparse.FileType('r'), required=True)
	parser.add_argument('-t','--tanimoto', help='Tanimoto similarity, 60 by default', type=int, default=60)
	parser.add_argument('-o','--output', help='Output result in a tanimoto.txt file', default=False)
	args = parser.parse_args()

	res = chembl_similarity(args.input, args.tanimoto)
	if not res:
		print('No result', file=sys.stderr)
		return 0

	df = pd.DataFrame.from_dict(json_normalize(res), orient='columns')
	df_final = df[['molecule_chembl_id', 'similarity','molecule_structures.canonical_smiles']]
	df_nodup= df_final.drop_duplicates(subset=['molecule_chembl_id'])
	output = StringIO()
	df_nodup.to_csv(output) ;

	if args.output:
		f = open(args.output, "w")
		f.write(str(df_nodup))
	else: print(output.getvalue())


if __name__=='__main__':
		Main()
