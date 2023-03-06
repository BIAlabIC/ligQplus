import argparse
from chembl_webresource_client.new_client import new_client
from io import StringIO
import pandas as pd
import numpy as np
import collections
from pandas import json_normalize
import sys

def chembl_trg (handle):
    trg_info=[]
    if handle:
        for molecule_chembl_id in handle:
            molecule_chembl_id = molecule_chembl_id.strip()
            trg_1=list(new_client.activity.filter(molecule_chembl_id=molecule_chembl_id,
            assay_type = 'B').only('target_chembl_id', 'activity_comment', 'pchembl_value', 'molecule_chembl_id'))
            trg_2=list(new_client.mechanism.filter(molecule_chembl_id =trg_id, max_phase__in = [3, 4]
            ).only("target_chembl_id", "max_phase"))
            trg_info=trg_info+trg_1+trg_2
    return trg_info

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Input tanimoto_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o','--output', help='Output result in a target.txt file',
     type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    records=chembl_trg(args.input)
    if not records:
    	print('No result', file=sys.stderr)
    	return 0

    df = pd.DataFrame.from_dict(json_normalize(records), orient='columns')
    df.loc[((df['activity_comment'] == 'Active') & (df['pchembl_value'].isnull())), 'pchembl_value'] = 6  #TODO: revisar este filtro
    df.loc[((df['max_phase'].notnull())), 'pchembl_value'] = 6
    df1 = df.dropna(subset=['pchembl_value'])
    def pchembl_median(x):
        names = {
            'activity_comment': x.iloc[0]['activity_comment'],
            'max_phase': x.iloc[0]['max_phase'],
            'pchembl_median': x['pchembl_value'].median(),
        }
        return(pd.Series(names, index = ['pchembl_median', 'activity_comment', 'max_phase']))
    df_pchembl_median = df1.groupby(['target_chembl_id']).apply(pchembl_median).reset_index()
    df_drop = df_pchembl_median.drop(df_pchembl_median[df_pchembl_median.pchembl_median < 6].index)
    df_nodup= df_drop.drop_duplicates(subset=['target_chembl_id'])
    output=StringIO()
    df_nodup.to_csv(args.output, sep='\t') ;


if __name__=='__main__':
		Main()
