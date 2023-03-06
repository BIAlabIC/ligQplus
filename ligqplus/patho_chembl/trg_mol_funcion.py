from chembl_webresource_client.new_client import new_client
import pandas as pd
from pandas import json_normalize
import json
import argparse
from io import StringIO
import sys

def trg_mol(chembl_trg_id):
    chembl_trg_id = chembl_trg_id.strip()
    mol_1=list(new_client.activity.filter(target_chembl_id =chembl_trg_id,
     assay_type = 'B').only("target_chembl_id", "molecule_chembl_id", "canonical_smiles", "pchembl_value",
     "target_organism", "activity_comment"))
    mol_2=list(new_client.mechanism.filter(target_chembl_id =chembl_trg_id, max_phase__in = [3, 4]
    ).only("molecule_chembl_id", "max_phase"))
    return mol_1, mol_2

def trg_mol_list (handle):
    mol_=[]
    if handle:
        for trg_id in handle:
            mol_1, mol_2 = trg_mol(trg_id)
            mol_=mol_+mol_1+mol_2
    return mol_

def pros_trg_mol_list(records):
     df = pd.DataFrame.from_dict(json_normalize(records), orient='columns')

     if 'max_phase' not in df.columns:
         df['max_phase'] = None
     elif 'pchembl_value' not in df.columns:
         df['pchembl_value'] = None
     elif 'activity_comment' not in df.columns:
         df['activity_comment'] = None
     df.loc[((df['activity_comment'] == 'Active') & (df['pchembl_value'].isnull())), 'pchembl_value'] = 6  #TODO: revisar este filtro |
     df.loc[((df['max_phase'].notnull())), 'pchembl_value'] = 6
     df.loc[((df['pchembl_value'].isnull())), 'pchembl_value'] = 0
     df1 = df.dropna(subset=['pchembl_value'])
     def pchembl_median(x):
         names = {
             'activity_comment': x.iloc[0]['activity_comment'],
             'max_phase': x.iloc[0]['max_phase'],
             'target_organism': x.iloc[0]['target_organism'],
             'pchembl_median': x['pchembl_value'].median()
         }
         return(pd.Series(names, index = ['pchembl_median', 'activity_comment', 'target_organism', 'max_phase']))
     df_pchembl_median = df1.groupby(['molecule_chembl_id']).apply(pchembl_median).reset_index()
     df_drop = df_pchembl_median.drop(df_pchembl_median[df_pchembl_median.pchembl_median < 6].index)
     df_nodup= df_drop.drop_duplicates(subset=['molecule_chembl_id'])
     return df_nodup

def trg_to_mol(chembl_trg_id):
    assay, mechanism = trg_mol(chembl_trg_id)
    records= assay+mechanism
    if records:
        return pros_trg_mol_list(records)
    else:
        return pd.DataFrame()

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Input. File containing list of target chembl_id', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o','--output', help='Output result must be a .csv file',
     type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    records=trg_mol_list(args.input)
    if records:
        df_nodup = pros_trg_mol_list(records)
        output=StringIO()
        df_nodup.to_csv(args.output, sep='\t') ;

    else:
    	print('No result', file=sys.stderr)
    	return 0

if __name__=='__main__':
		Main()
