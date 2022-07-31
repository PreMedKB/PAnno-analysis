from panno import genotype_resolution
import pandas as pd
import numpy as np
import os, re

gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",\
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", \
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
########### Test the accuracy
pharmgkb_100genome = {'AFR': 'African American/Afro-Caribbean', 'AMR': 'Latino', 'EAS': 'East Asian', 'EUR': 'European'}
metadata = pd.read_csv('../res/population.txt', sep='\t')
test_vcfs = os.listdir('../vcf/')
extract_match_res = {}
final_rank_res = {}
n = 0
for vcf in test_vcfs:
  if vcf.endswith('.vcf'):
    n = n+1
    print("%i/88" % n)
    vcf_file = '../vcf/%s' % vcf
    # Get sample and race info
    sample = vcf.split('.')[0]
    race = metadata.loc[metadata['Get-RM 137 Samples']==sample, 'Superpopulation'].to_list()[0]
    race = pharmgkb_100genome[race]
    print(sample, race)
    res = genotype_resolution.predict_diplotype.main(vcf_file, race, gene_list)
    # Extract match res
    extract_match_res[sample] = {}
    for gene in gene_list:
      extract_match_res[sample][gene] = res[gene]['extract_match_res']
    # Final rank res
    final_rank_res[sample] = {}
    for gene in gene_list:
      final_rank_res[sample][gene] = res[gene]['final_rank_res']


extract_df = pd.DataFrame(extract_match_res)
final_df = pd.DataFrame(final_rank_res)
extract_df.T.to_csv('../manuscript/diplotype_concordance/detail_res1_v20220420.txt', sep='\t')
final_df.T.to_csv('../manuscript/diplotype_concordance/detail_res2_v20220420.txt', sep='\t')



### Validation data
getrm = pd.read_csv('../cpat_data/test/getrm_raw.txt', sep='\t', index_col=0)
# Filter samples and genes
samples = list(final_rank_res.keys())
overlap_gene = list(set(gene_list).intersection(set(getrm.columns.to_list())))
getrm_new = getrm[getrm.index.isin(samples)][overlap_gene].T.replace(np.nan, '-')

pharmcat = pd.read_csv('../cpat_data/test/pharmcat_raw.txt', sep='\t', index_col=0)
# Filter samples and genes
samples = sorted(list(final_rank_res.keys()))
pharmcat_new = pharmcat[pharmcat.index.isin(samples)][overlap_gene].T.replace(np.nan, '-')
for s in ['NA10856', 'NA10846', 'NA12336', 'NA12753', 'NA10855', 'NA12236', 'NA07019', 'NA19109', 'NA18484', 'NA07055', 'NA10854', 'NA18518', 'NA11993', 'NA18855', 'NA07348', 'NA11839', 'NA19176', 'NA19226', 'NA19122', 'NA06993', 'NA12145', 'NA10859', 'NA06991', 'NA19174', 'NA10865', 'NA10838', 'NA07029', 'NA12892', 'NA10831']:
  pharmcat_new[s] = np.nan

# # samples used in PharmCAT
# pharmcat_samples = []
# for i in metadata['PharmCAT 59 Samples'].to_list():
#   if i != '-' and pd.isna(i) is False:
#     pharmcat_samples.append(i)

# # Accuracy by gene
# for g in overlap_gene:
#   test = test_df.loc[g, samples].to_list()
#   validate = getrm_new.loc[g, samples].to_list()
#   # remove the nan
#   index = [x for x, y in list(enumerate(validate)) if y != '-']
#   right = 0
#   for i in index:
#     if g in ['ABCG2', 'VKORC1']:
#       x = test[i]
#       genotype = validate[i].replace('G', 'C').replace('A', 'T')
#       yy = [genotype[0], genotype[1]]
#     else:
#       if g == 'DPYD':
#         x = (test[i][0].replace('Reference', '*1'), test[i][1].replace('Reference', '*1'))
#       else:
#         x = test[i]
#       y = re.findall(r'\*\w+', validate[i])
#       yy = [j.replace('*0', '*') for j in y]
    
#     if set(x) <= set(yy):
#       right = right+1
#     else:
#       # Has the overlap of haplotype is passed
#       if 'no consensus' in validate[i]:
#         y = re.findall(r'\*\w+', validate[i])
#         yy = [j.replace('*0', '*') for j in y]
#         if len(set(x).intersection(set(yy))) > 0:
#           right = right+1
  
#   print(g, round(right/len(index)*100, 2))



###################################
###### Supplementary table 1 ######
###################################
output = []
index_name = []
final_df = pd.DataFrame(final_rank_res)
for g in sorted(overlap_gene):
  test = final_df.loc[g, samples].to_list()
  validate = getrm_new.loc[g, samples].to_list()
  pharmcat = pharmcat_new.loc[g, samples].to_list()
  cpat_new = test
  # for item in test:
  #   item_new = []
  #   for tmp in item:
  #     i = tmp.split('/')
  #     ii = i[0].replace('Reference', '*1') + '/' + i[1].replace('Reference', '*1')
  #     if g in ['ABCG2', 'VKORC1', 'IFNL3']:
  #       ii = i[0].replace('C', 'G').replace('T', 'A') + i[1].replace('C', 'G').replace('T', 'A')
  #     item_new.append(ii)
  #   cpat_new.append(';'.join(item_new))
  validate_new = [i.replace(' ', '') for i in validate]
  index_name.extend(['%s.GetRM' % g, '%s.PharmCAT' % g, '%s.CPAT' % g])
  output.extend([validate_new, pharmcat, cpat_new])

output_df = pd.DataFrame(output, index = index_name, columns = samples).T
output_df.to_csv('../manuscript/diplotype_concordance/GetRM_CPAT_PharmCAT_v20220421.txt', sep = '\t')



###################################
###### Supplementary table 2 ######
###################################
loci = pd.read_csv('../Get-RM/loci_alleles_tested.txt', sep='\t')
lt = list()
for index, row in loci.iterrows():
  gene = row[0]
  assays = list(row[1:].dropna().values)
  assays.append('*1')
  uniq_loci = set(', '.join(assays).split(', '))
  # cpat
  allele_definition_table = "./assets/definition/%s_allele_definition_table.txt" % (gene)
  cpat_gene_define = pd.read_csv(allele_definition_table, sep='\t', skiprows=5)
  print(cpat_gene_define.columns[0])
  cpat_loci = set(cpat_gene_define[cpat_gene_define.columns[0]].values)
  getrm_uniq = uniq_loci - cpat_loci
  cpat_uniq = cpat_loci - uniq_loci
  lt.append([gene, ', '.join(sorted(uniq_loci)), ', '.join(sorted(cpat_loci)), ', '.join(sorted(getrm_uniq)), ', '.join(sorted(cpat_uniq))])

df = pd.DataFrame(lt)
df.to_csv('../Get-RM/loci_alleles_uniq.txt', sep='\t')