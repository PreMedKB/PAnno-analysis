#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Data list for dosing commendations.
"""

import sqlite3
import pandas as pd
import numpy as np

## Connected database
pgx_kb_fp = '../PAnno-package/panno/assets/pgx_kb.sqlite3'
conn = sqlite3.connect(pgx_kb_fp)
cursor = conn.cursor()

dosing = cursor.execute('SELECT DISTINCT Source, RelatedGeneID, RelatedDrugID FROM ClinDosingGuideline;')
dosing_df = pd.DataFrame(dosing.fetchall(), columns=['Source', 'RelatedGeneID', 'RelatedDrugID'])
anno = cursor.execute('SELECT DISTINCT Gene, Drug, GeneID, DrugID FROM ClinAnn WHERE EvidenceLevel != 3 AND EvidenceLevel != 4;')
anno_df = pd.DataFrame(anno.fetchall(), columns=['Gene', 'Drug', 'RelatedGeneID', 'RelatedDrugID'])

filter_dosing_df = pd.merge(dosing_df, anno_df).drop_duplicates()