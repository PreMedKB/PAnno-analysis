from panno import genotype_resolution, clinical_annotation, pgx_report
import pandas as pd
import os

gene_list = ["G6PD", "MT-RNR1", "ABCG2", "CACNA1S", "CFTR", "IFNL3", "VKORC1", "RYR1", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "NUDT15", "SLCO1B1", "TPMT", "UGT1A1"]


########### Test the accuracy
metadata = pd.read_csv('./metadata.txt', sep='\t')
test_vcfs = os.listdir('./vcf/')
outdir = './report/'
# Diplotypes
exact_res = {}; step1_res = {}; step2_res = {}
# Annotations
summary = {}
prescription_df = pd.DataFrame()
phenotype_df = pd.DataFrame()
clinical_anno_df = pd.DataFrame()

for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    germline_vcf = './vcf/%s' % vcf
    # Get sample and race info
    sample_id = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==sample_id, 'Abbreviation'].to_list()[0]
    population = race.upper()
    pop_dic = {'AAC': 'African American/Afro-Caribbean', 'AME': 'American', 'SAS': 'Central/South Asian', 'EAS': 'East Asian', 'EUR': 'European', 'LAT': 'Latino', 'NEA': 'Near Eastern', 'OCE': 'Oceanian', 'SSA': 'Sub-Saharan African'}
    if population not in pop_dic.keys():
      print('The input population is not included in PAnno. Please check if the abbreviation is used correctly.')
    race = pop_dic[race]
    print(sample_id, race)
    ### Run PAnno
    print('  - Parsing PGx related genotypes ...')
    dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(race, germline_vcf)
    ##### Diplotypes concordance
    # Exact match res
    exact_res[sample_id] = {}
    for gene in gene_list:
      exact_res[sample_id][gene] = dic_diplotype[gene]['exact_res']
    # Step1 rank res
    step1_res[sample_id] = {}
    for gene in gene_list:
      step1_res[sample_id][gene] = dic_diplotype[gene]['step1_res']
    # Final rank res
    step2_res[sample_id] = {}
    for gene in gene_list:
      step2_res[sample_id][gene] = dic_diplotype[gene]['step2_res']
    
    print('  - Annotating clinical information ...')
    summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    print('  - Generating PAnno report ...')
    race = "%s (%s)" % (pop_dic[population], population)
    fp = os.path.join(outdir, "%s.PAnno.html" % sample_id)
    pgx_report.report(race, summary, prescribing_info, multi_var, single_var, phenotype_predict, clinical_anno, fp, sample_id)
    
    ##### Clinical annotation results
    summary[sample_id] = summary
    prescribing_info['Sample'] = sample_id
    phenotype_predict['Sample'] = sample_id
    clinical_anno['Sample'] = sample_id
    prescription_df = pd.concat([prescription_df, prescribing_info], axis=0)
    phenotype_df = pd.concat([phenotype_df, phenotype_predict], axis=0)
    clinical_anno_df = pd.concat([clinical_anno_df, clinical_anno], axis=0)


# Get all drugs
import sqlite3
pgx_kb_fp = "../PAnno/panno/assets/pgx_kb.sqlite3"
conn = sqlite3.connect(pgx_kb_fp)
cursor = conn.cursor()
guide = cursor.execute("SELECT Drug FROM GuidelineMerge WHERE Source NOT IN ('AusNZ', 'SEFF', 'ACR', 'CFF');")
guide = cursor.fetchall()
summary_drug = []
for drug in guide:
  summary_drug.append(drug[0])


summary_dic = {}
for s in summary.keys():
  single = {}
  for d in summary_drug:
    if d in summary[s]['Avoid']:
      single[d] = 'Avoid'
    elif d in summary[s]['Caution']:
      single[d] = 'Caution'
    elif d in summary[s]['Routine']:
      single[d] = 'Routine'
    else:
      single[d] = 'Missing'
  summary_dic[s] = single

summary_df = pd.DataFrame(summary_dic)
summary_df.to_csv('./res/summary_df.txt', sep="\t")
prescription_df.to_csv('./res/prescription_df.txt', sep="\t", index=0)
phenotype_df.to_csv('./res/phenotype_df.txt', sep="\t", index=0)
clinical_anno_df.to_csv('./res/clinical_anno_df.txt', sep="\t", index=0)

exact_df = pd.DataFrame(exact_res).T
step1_df = pd.DataFrame(step1_res).T
step2_df = pd.DataFrame(step2_res).T

## Merge the above dfs into one df
exact_df.columns = ['%s.exact' % i for i in exact_df.columns]
step1_df.columns = ['%s.step1' % i for i in step1_df.columns]
step2_df.columns = ['%s.step2' % i for i in step2_df.columns]
df = pd.concat([exact_df, step1_df, step2_df], axis=1)
df.loc[sorted(df.index), sorted(df.columns)].to_csv('./res/diplotype_different_steps.txt', sep='\t')




ann = cursor.execute("SELECT * FROM ClinAnn WHERE EvidenceLevel != '3';")
ann = cursor.fetchall()
ann_df = pd.DataFrame(ann, columns=['ID', 'CAID', 'Gene', 'Variant', 'Allele1', 'Allele2', 'Annotation1', 'Annotation2', 'Function1', 'Function2', 'Score1', 'Score2', 'CPICPhenotype', 'PAnnoPhenotype', 'Drug', 'Phenotypes', 'EvidenceLevel', 'LevelOverride', 'LevelModifier', 'Score', 'PMIDCount', 'EvidenceCount', 'Specialty', 'PhenotypeCategory'])
ann_df.insert(ann_df.shape[1], 'Allele', '')

for index, row in ann_df.iterrows():
  # Update Variant
  if ', ' in row.Variant:
    # variants = row.Variant.split(', %s' % row.Gene)
    # for i in range(1, len(variants)):
    #     variants[i] = row.Gene + variants[i]
    # row.Variant = '%'.join(variants)
    # row.Variant = row.Gene+'**'
    row.Variant = row.Gene+row.Allele1+' % '+row.Gene+row.Allele2
  # Add Allele
  if row.Variant.startswith('rs'):
    row.Allele = row.Gene + row.Variant + row.Allele1 + row.Allele2
  else:
    row.Allele = row.Gene + row.Allele1 + row.Allele2
  # Update this row
  ann_df.iloc[index] = row


ann_df = ann_df.drop(['Variant'], axis=1).join(ann_df['Variant'].str.split(' % ', expand=True).stack().reset_index(level=1, drop=True).rename('Variant')).drop_duplicates().reset_index(drop = True)
ann_df.to_csv('./res/raw_anno_v20221212.txt', sep="\t", index=0)