from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from importlib import reload
import collections
from pandas import json_normalize
import json
import argparse
import sys
from sqlalchemy import create_engine
import sqlite3
from importlib import reload
import os

def search_bypfam(dataset_path):
    engine = create_engine(dataset_path)
    CHEMBL_VERSION = 27
    find_molbypfam = ('''SELECT td.chembl_id as target_chemblid, a2.pchembl_value,
    a2.activity_comment, md.chembl_id as compound_chemblid, source_domain_id,
    dm.mec_id, max_phase
    FROM drug_mechanism dm
    JOIN binding_sites bs on bs.tid = dm.tid
    JOIN target_dictionary td ON td.tid = bs.tid
    JOIN site_components sc ON sc.site_id =bs.site_id
    JOIN domains d2 ON d2.domain_id = sc.domain_id
    JOIN activities a2 ON dm.molregno = a2.molregno
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN compound_properties cp ON cp.molregno = md.molregno
    JOIN compound_records cr ON cr.molregno = cp.molregno
    WHERE dm.mec_id IS NOT NULL
    AND a2.activity_comment LIKE 'Active'
    UNION
    SELECT td.chembl_id as target_chemblid, a2.pchembl_value, a2.activity_comment,
    md.chembl_id as compound_chemblid, source_domain_id, dm.mec_id, max_phase
    FROM drug_mechanism dm
    JOIN binding_sites bs on bs.tid = dm.tid
    JOIN target_dictionary td ON td.tid = bs.tid
    JOIN site_components sc ON sc.site_id =bs.site_id
    JOIN domains d2 ON d2.domain_id = sc.domain_id
    JOIN activities a2 ON dm.molregno = a2.molregno
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN compound_properties cp ON cp.molregno = md.molregno
    JOIN compound_records cr ON cr.molregno = cp.molregno
    WHERE a2.src_id = 15
    AND dm.mec_id IS NOT NULL
    AND a2.pchembl_value >= 6.0
    AND cp.PSA IS NOT NULL
    UNION
    SELECT td.chembl_id as target_chemblid, a2.pchembl_value, a2.activity_comment,
    md.chembl_id as compound_chemblid, source_domain_id, dm.mec_id, max_phase
    FROM drug_mechanism dm
    JOIN binding_sites bs on bs.tid = dm.tid
    JOIN target_dictionary td ON td.tid = bs.tid
    JOIN site_components sc ON sc.site_id =bs.site_id
    JOIN domains d2 ON d2.domain_id = sc.domain_id
    JOIN activities a2 ON dm.molregno = a2.molregno
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN compound_properties cp ON cp.molregno = md.molregno
    JOIN compound_records cr ON cr.molregno = cp.molregno
    WHERE md.max_phase >= 3.0
    AND dm.mec_id IS NOT NULL; ''')

    df_mol = pd.read_sql(find_molbypfam, engine)
    return df_mol


def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--dataset", help="Path to the directory of ChEMBL DB",
    default= "chembl_27.db")
    parser.add_argument('-o','--output', help='Output result must be .csv file',
     type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()

    db='sqlite:///'+ os.path.abspath(args.dataset)

    if args.dataset:
        print('Database is being generated. This may take 5-10min', file=sys.stderr)
        df_targets=search_bypfam(db)
        df_targets.to_csv("/tmp/mech.csv")
        if len(df_targets):
            df_drop=df_targets.drop_duplicates(subset= ['target_chemblid', 'source_domain_id', 'compound_chemblid'])
            grouped_df = df_drop.groupby(["target_chemblid", "compound_chemblid"])
            grouped_pfam = grouped_df.agg({'source_domain_id':list, 'max_phase':'max'})
            grouped_pfam = grouped_pfam.reset_index()
            grouped_pfam['Domain_key']=['|'.join(sorted(set(x))) for x in grouped_pfam.source_domain_id]
            result = grouped_pfam
            result.to_csv(args.output) ;
    else:
        print('No database', file=sys.stderr)

    return 0

if __name__=='__main__':
		Main()
