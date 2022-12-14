import pandas as pd
import numpy as np
import re, itertools, json

panno_dip_base = json.loads(open("./assets/panno_dip_base.json").read())

gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",\
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", \
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]

for gene in gene_list:
  info = panno_dip_base[gene]
  ## Combine the haplotypes to get diplotypes
  hap_define = info['haplotype_definition']
  
dic_degenerate_bases = {
  'R': ['A', 'G'],
  'Y': ['C', 'T'],
  'M': ['A', 'C'],
  'K': ['G', 'T'],
  'S': ['G', 'C'],
  'W': ['A', 'T'],
  'H': ['A', 'T', 'C'],
  'B': ['G', 'T', 'C'],
  'V': ['G', 'A', 'C'],
  'D': ['G', 'A', 'T'],
  'N': ['A', 'T', 'C', 'G']
}

dic_gene_ref_haplotype = {
  'ABCG2': 'G',
  'CACNA1S': 'Reference',
  'CYP2B6': '*1',
  'CYP2C19': '*38',
  'CYP2C9': '*1',
  'CYP2C8': '*1',
  'CYP2D6': '*1',
  'CYP3A4': '*1',
  'CYP3A5': '*1',
  'CYP4F2': '*1',
  'CFTR': 'ivacaftor non-responsive CFTR sequence',
  'DPYD': 'Reference',
  'MT-RNR1': 'Reference',
  'G6PD': 'B (wildtype)',
  'NUDT15': '*1',
  'RYR1': 'Reference',
  'SLCO1B1': '*1',
  'TPMT': '*1',
  'UGT1A1': '*1',
  'VKORC1': 'C',
  'IFNL3': 'C'
}

statistic = []
for gene in gene_list:
  allele_definition_table = "./assets/definition/%s_allele_definition_table.txt" % (gene)
  pos = []
  list_alt = []
  list_ref = []
  index_nonsnp = []
  dic_alle2genotype = {}; dic_alle2base = {}
  with open(allele_definition_table, 'r', encoding='utf-8') as file:
    flag = 0
    for line in file:
      line = line.replace("\n", "")
      info = line.split('\t')
      flag_judge = info[0]
      if 'GRCh38' in line:
        for i in info:
          matchobj = re.search(r'\w\.(\d+)(\w)\>(\w)', i)
          #matchobj = re.search(r'g\.(\d+)(\w)\>(\w)', i)
          if matchobj:
            pos.append(matchobj.group(1))
            list_alt.append(matchobj.group(2))
          else:
            matchobj = re.search(r'\w\.(\d+)(\_\d+)?(del|ins)(\w*)', i)
            #matchobj = re.search(r'g\.(\d+)del(\w*)', i)
            if matchobj:
              pos.append(matchobj.group(1))
              list_alt.append(matchobj.group(4))
            elif "GRCh38" not in i:
              index_nonsnp.append(info.index(i)-1)
              pos.append(info.index(i)-1)
              list_alt.append(info.index(i)-1)
      
      else:
        if "rsID" in line:
          for i in index_nonsnp:
            pos[i] = info[i+1]
            list_alt[i] = info[i+1]
      
      # Start to process the haplotypes
      if flag == 1:
        alle = info.pop(0)
        
        if alle == dic_gene_ref_haplotype[gene]:
          list_ref = info
          gt = [0 for genotype in range(len(info))]
          define = info
        else:
          gt = []
          define = []
          for j in range(0, len(info)):
            if info[j] == '' or info[j] == list_ref[j]:
              genotype = 0
              allele = list_ref[j]
            elif info[j] in list(dic_degenerate_bases.keys()):
              genotype = []
              allele = dic_degenerate_bases[info[j]]
              for wobble in dic_degenerate_bases[info[j]]:
                if wobble == list_ref[j]:
                  genotype.append(0)
                else:
                  genotype.append(1)
            else:
              allele = info[j]
              genotype = 1
            gt.append(genotype)
            define.append(allele.split('; '))
        dic_alle2genotype[alle] = gt
        dic_alle2base[alle] = define
      
      if re.search('%s Allele' % gene, flag_judge, re.IGNORECASE):
        flag = 1
  
  ## 基于dic_alle2base进行两两组合
  import itertools
  # candidate = {}
  combination = list(itertools.combinations_with_replacement(dic_alle2base, 2))
  print(gene, len(combination))
  loci_count = len(list(dic_alle2base.values())[0])
  diplotype_df = pd.DataFrame()
  for comb in combination:
    diplotype_define = []
    i = comb[0]; j = comb[1]
    for index in range(0, loci_count):
      hap1 = dic_alle2base[i][index]
      hap2 = dic_alle2base[j][index]
      if type(hap1) is str:
        hap1 = [hap1]
      if type(hap2) is str:
        hap2 = [hap2]
      # cut_gt_tmp
      cut_gt_tmp = []
      for item in itertools.product(hap1, hap2):
        cut_gt_tmp.append("/".join(sorted(item)))
      diplotype_define.append("; ".join(cut_gt_tmp))
    ## Split by ;
    df = pd.DataFrame(diplotype_define)
    df = df[0].str.split('; ', expand=True)
    for index, row in df.iterrows():
      for k in range(0, len(row)):
        if row[k] == None:
          row[k] = row[0]
      df.iloc[index, ] = row
    df.columns = [i+'/'+j]*len(row)
    diplotype_df = pd.concat([diplotype_df, df], axis=1)
  
  test_df = diplotype_df.T
  # Gene, haplotype length, diplotype length, duplicated diplotypes
  statistic.append([gene, len(dic_alle2base.keys()), test_df.shape[0], test_df.duplicated(keep=False).sum()])
  print(statistic)
  test_df['Duplicated'] = test_df.duplicated(keep=False)
  test_df_duplicate = test_df[test_df['Duplicated'] == True]
  test_df_duplicate = test_df_duplicate.sort_values(test_df_duplicate.columns.to_list())
  test_df_duplicate.to_csv('./same_pattern/%s_diplotype.txt' % gene, sep="\t")

same_pattern = pd.DataFrame(statistic, columns=['Gene', 'Haplotype', 'Diplotype', 'DuplicatedDiplotype'])
same_pattern.to_csv('same_pattern.txt', sep='\t', index=0)