from panno import genotype_resolution, clinical_annotation, pgx_report
import pandas as pd
import os

gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",\
            "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", \
            "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]

########### Test the accuracy
# race_dic = {'AFR': 'African American/Afro-Caribbean', 'AMR': 'Latino', 'EAS': 'East Asian', 'EUR': 'European'}
metadata = pd.read_csv('../metadata.txt', sep='\t')
test_vcfs = os.listdir('../vcf/')
outdir = '../report/'

clinical_anno_df = pd.DataFrame()
dosing_guideline_df = pd.DataFrame()
pgx_summary_df = pd.DataFrame()
exact_res = {}
step1_res = {}
step2_res = {}

for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    germline_vcf = './data/vcf/%s' % vcf
    # Get sample and race info
    sample_id = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==sample_id, 'Abbreviation'].to_list()[0]
    population = race.upper()
    pop_dic = {'AAC': 'African American/Afro-Caribbean', 'AME': 'American', 'SAS': 'Central/South Asian', 'EAS': 'East Asian', 'EUR': 'European', 'LAT': 'Latino', 'NEA': 'Near Eastern', 'OCE': 'Oceanian', 'SSA': 'Sub-Saharan African'}
    if population not in pop_dic.keys():
      print('The input population is not included in CPAT. Please check if the abbreviation is used correctly.')
    race = pop_dic[race]
    print(sample_id, race)
    ### Run CPAT
    print('  - Parsing PGx related genotypes ...')
    dic_diplotype, dic_rs2gt, hla_subtypes = genotype_resolution.resolution(race, germline_vcf)
    print('  - Annotating clinical information ...')
    pgx_summary, clinical_anno_table, dosing_guideline_table = clinical_annotation.annotation(dic_diplotype, dic_rs2gt, hla_subtypes)
    print('  - Generating CPAT report ...')
    pgx_report.report(race, pgx_summary, dic_diplotype, clinical_anno_table, dosing_guideline_table, outdir, sample_id)
    ### Collect the results
    # Clinical annotation
    pgx_summary['Sample'] = sample_id; pgx_summary['Drug'] = pgx_summary.index
    clinical_anno_table['Sample'] = sample_id
    dosing_guideline_table['Sample'] = sample_id
    pgx_summary_df = pd.concat([pgx_summary_df, pgx_summary], axis=0)
    clinical_anno_df = pd.concat([clinical_anno_df, clinical_anno_table], axis=0)
    dosing_guideline_df = pd.concat([dosing_guideline_df, dosing_guideline_table], axis=0)
    # Diplotypes
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


clinical_anno_df.EvidenceLevel.drop_duplicates()
pgx_summary_df.groupby('EvidenceLevel').count()

clinical_anno_df.to_csv('../res/clinical_anno_df.txt', sep="\t", index=0)
dosing_guideline_df.to_csv('../res/dosing_guideline_df.txt', sep="\t", index=0)
pgx_summary_df.to_csv('../res/pgx_summary_df.txt', sep="\t", index=0)

exact_df = pd.DataFrame(exact_res).T
step1_df = pd.DataFrame(step1_res).T
step2_df = pd.DataFrame(step2_res).T

## Merge the above dfs into one df
exact_df.columns = ['%s.exact' % i for i in exact_df.columns]
step1_df.columns = ['%s.step1' % i for i in step1_df.columns]
step2_df.columns = ['%s.step2' % i for i in step2_df.columns]
df = pd.concat([exact_df, step1_df, step2_df], axis=1)
df[sorted(df.columns)].to_csv('../res/diplotype_different_steps.txt', sep='\t')
