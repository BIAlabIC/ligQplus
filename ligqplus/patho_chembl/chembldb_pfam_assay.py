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

#def load_dataset(dataset_path):
#    chembl = dataset_path
#    engine = create_engine(dataset_path)
#    CHEMBL_VERSION = 27
#    return engine

def search_bypfam(dataset_path):
    engine = create_engine(dataset_path)
    CHEMBL_VERSION = 27
    find_molbypfam = ('''SELECT g.chembl_id as target_chemblid, source_domain_id,
    c.chembl_id as compound_chemblid, a.pchembl_value, a.activity_comment
    FROM  ACTIVITIES a
    JOIN MOLECULE_HIERARCHY b ON a.MOLREGNO = b.MOLREGNO
    JOIN MOLECULE_DICTIONARY c ON b.PARENT_MOLREGNO = c.MOLREGNO
    JOIN COMPOUND_PROPERTIES d ON b.PARENT_MOLREGNO = d.MOLREGNO
    JOIN COMPOUND_STRUCTURES e ON b.PARENT_MOLREGNO = e.MOLREGNO
    JOIN ASSAYS f ON a.ASSAY_ID = f.ASSAY_ID
    JOIN TARGET_DICTIONARY g ON f.TID = g.TID
    JOIN DOCS h ON a.DOC_ID = h.DOC_ID
    JOIN TARGET_COMPONENTS i ON g.tid = i.tid
    JOIN COMPONENT_CLASS j ON i.component_id = j.component_id
    JOIN PROTEIN_CLASSIFICATION k ON j.protein_class_id = k.protein_class_id
    JOIN component_domains ON component_domains.component_ID = i.component_ID
    JOIN domains ON domains.domain_ID = component_domains.domain_ID
    JOIN compound_records ON a.record_ID = compound_records.record_ID
    WHERE a.pchembl_value >= 6.0
    AND (a.potential_duplicate = 0 OR a.potential_duplicate IS NULL)
    AND assay_type = 'B'
    UNION
    SELECT g.chembl_id as target_chemblid, source_domain_id, c.chembl_id as compound_chemblid,
    a.pchembl_value, a.activity_comment
    FROM  ACTIVITIES a
    JOIN MOLECULE_HIERARCHY b ON a.MOLREGNO = b.MOLREGNO
    JOIN MOLECULE_DICTIONARY c ON b.PARENT_MOLREGNO = c.MOLREGNO
    JOIN COMPOUND_PROPERTIES d ON b.PARENT_MOLREGNO = d.MOLREGNO
    JOIN COMPOUND_STRUCTURES e ON b.PARENT_MOLREGNO = e.MOLREGNO
    JOIN ASSAYS f ON a.ASSAY_ID = f.ASSAY_ID
    JOIN TARGET_DICTIONARY g ON f.TID = g.TID
    JOIN DOCS h ON a.DOC_ID = h.DOC_ID
    JOIN TARGET_COMPONENTS i ON g.tid = i.tid
    JOIN COMPONENT_CLASS j ON i.component_id = j.component_id
    JOIN PROTEIN_CLASSIFICATION k ON j.protein_class_id = k.protein_class_id
    JOIN component_domains ON component_domains.component_ID = i.component_ID
    JOIN domains ON domains.domain_ID = component_domains.domain_ID
    JOIN compound_records ON a.record_ID = compound_records.record_ID
    WHERE  a.activity_comment LIKE 'Active'
    AND (a.potential_duplicate = 0 OR a.potential_duplicate IS NULL)
    AND assay_type = 'B'
    AND g.target_type IN ('SINGLE PROTEIN', 'PROTEIN COMPLEX');  ''')
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
        if len(df_targets):
            df_drop=df_targets.drop_duplicates(subset= ['target_chemblid', 'source_domain_id', 'compound_chemblid'])
            grouped_df = df_drop.groupby("target_chemblid")
            grouped_pfam = grouped_df["source_domain_id"].apply(list)
            grouped_df = df_drop.groupby("target_chemblid")
            grouped_pfam = grouped_df["source_domain_id"].apply(list)
            grouped_pfam = grouped_pfam.reset_index()
            grouped_pfam['Domain_key']=['_'.join(sorted(set(x))) for x in grouped_pfam.source_domain_id]
            grouped_compound = grouped_df["compound_chemblid"].apply(list)
            grouped_compound = grouped_compound.reset_index()
            result = pd.merge(grouped_compound, grouped_pfam,  how='left', on=['target_chemblid'])
            result.to_csv(args.output) ;
    else:
        print(f'No database', file=sys.stderr)

    return 0

#vscodesqlite3.connect

if __name__=='__main__':
		Main()
