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
    find_pfam = (f''' SELECT a2.pchembl_value, a2.activity_comment, md.chembl_id as mol_chemblid,
    cr.compound_name, source_domain_id,  md.max_phase
    FROM drug_mechanism dm
    JOIN binding_sites bs on bs.tid = dm.tid
    JOIN site_components sc ON sc.site_id =bs.site_id
    JOIN domains d2 ON d2.domain_id = sc.domain_id
    JOIN activities a2 ON dm.molregno = a2.molregno
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    JOIN compound_properties cp ON cp.molregno = md.molregno
    JOIN compound_records cr ON cr.molregno = cp.molregno
    WHERE a2.src_id = 15
    AND a2.standard_type = 'IC50'
    AND source_domain_id ="{pfam_id}"
    AND (NOT (a2.activity_comment LIKE 'inconclusive' OR a2.activity_comment LIKE 'undetermined')
    OR a2.activity_comment IS NULL)
    AND cp.PSA IS NOT NULL;''')
#    print(find_pfam)
    df = pd.read_sql(find_pfam, engine)
    return df
#se puede agregar , g.pref_name a la query en sql

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='Input pfam_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o','--output', help='Output result must be .csv file',
     type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    for pfam in args.input:
        pfam=pfam.strip()
        if pfam:
            df_drugs=search_bypfam(pfam)
            if len(df_drugs):
                df_drop = df_drugs.drop_duplicates("mol_chemblid")
                df_drop.to_csv(args.output) ;
            else:
                print(f'No result for {pfam}', file=sys.stderr)

    return 0



if __name__=='__main__':
		Main()
