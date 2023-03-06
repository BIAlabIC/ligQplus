import pandas as pd
from pandas import json_normalize
import json
import argparse
import sys
from sqlalchemy import create_engine
import sqlite3
from importlib import reload
chembl ="sqlite:////home/fleer/Desktop/Tesis/Chembldb/chembl_27/chembl_27_sqlite/chembl_27/chembl_27_sqlite/chembl_27.db"
engine = create_engine(chembl)
CHEMBL_VERSION = 27

def search_bypfam(pfam_id):
    find_molbypfam = (f''' SELECT a.pchembl_value, a.activity_comment, c.chembl_id as mol_chemblid,
    compound_records.compound_name, g.chembl_id as target_chemblid, source_domain_id, g.pref_name, c.max_phase
    FROM ACTIVITIES a
    JOIN MOLECULE_HIERARCHY b ON a.MOLREGNO = b.MOLREGNO
    JOIN MOLECULE_DICTIONARY c ON b.PARENT_MOLREGNO = c.MOLREGNO
    JOIN ASSAYS f ON a.ASSAY_ID = f.ASSAY_ID
    JOIN TARGET_DICTIONARY g ON f.TID = g.TID
    JOIN TARGET_COMPONENTS i ON g.tid = i.tid
    JOIN compound_properties cp ON cp.molregno = c.molregno
    JOIN component_domains ON component_domains.component_ID = i.component_ID
    JOIN domains ON domains.domain_ID = component_domains.domain_ID
    JOIN compound_records ON a.record_ID = compound_records.record_ID
    WHERE a.pchembl_value IS NOT NULL
    AND (a.potential_duplicate = 0 OR a.potential_duplicate IS NULL)
    AND a.data_validity_comment IS NULL
    AND ((f.confidence_score = 9 AND target_type = 'SINGLE PROTEIN')
    OR (f.confidence_score = 7 AND target_type = 'PROTEIN COMPLEX')
    OR (f.confidence_score = 8 AND target_type = 'SINGLE PROTEIN')
    OR (f.confidence_score = 6 AND target_type = 'PROTEIN COMPLEX'))
    AND f.src_id = 1
    AND (NOT (a.activity_comment LIKE 'inconclusive' OR a.activity_comment LIKE 'undetermined')
    OR a.activity_comment IS NULL)
    AND assay_type = 'B'
    AND source_domain_id = "{pfam_id}" AND cp.psa IS NOT NULL ''')
    df_mol = pd.read_sql(find_molbypfam, engine)
    return df_mol


def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Input pfam_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o','--output', help='Output result must be .csv file',
     type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    for pfam in args.input:
        pfam=pfam.strip()
        if pfam:
            df_targets=search_bypfam(pfam)
            if len(df_targets):
                df_drop = df_targets.drop_duplicates("mol_chemblid")
                df_drop.to_csv(args.output) ;
            else:
                print(f'No result for {pfam}', file=sys.stderr)

    return 0



if __name__=='__main__':
		Main()
