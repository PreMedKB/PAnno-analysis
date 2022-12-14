setwd('/Volumes/TM/PAnno/PAnno-analysis/')
library(ggplot2); library(ggsci); library(ggthemes); library(ggpubr); library(ComplexHeatmap); library(tidyverse)

summary_df <- read.csv('./res/summary_df.txt', sep="\t", row.names = 'X')
prescription_df <- read.csv('./res/prescription_df.txt', sep="\t")
phenotype_df <- read.csv('./res/phenotype_df.txt', sep="\t")
clinical_anno_df <- read.csv('./res/clinical_anno_df.txt', sep="\t")

#race_dic = {'AFR': 'African American/Afro-Caribbean', 'AMR': 'Latino', 'EAS': 'East Asian', 'EUR': 'European'}
meta <- read.csv('./res/population.txt', sep='\t')
meta <- meta[c('Get.RM.137.Samples', 'Group')]; colnames(meta) <- c('Sample', 'Group')

# The proportion of CYP genes in different populations is meaningless, and the sample size is too small
dip_df <- unique(prescription_df[,c('Sample', 'Gene', 'Phenotype')])
dip_df <- merge(meta, dip_df)

dip_df %>% filter(Gene %in% c('CYP2B6', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5')) %>%
ggplot(aes(x = Gene, fill = Phenotype)) + 
  #geom_bar(position = "fill") +
  geom_bar(position = position_dodge(preserve = "single"))+
  labs(y = "Proportion") + facet_grid(Group~.)


##### FIGURE: Ratio of drug/gene and drug/variant
raw_anno <- read.csv('./res/raw_anno_v20221212.txt', sep="\t") %>% 
  select(Gene, Drug, PhenotypeCategory, Variant, Allele1, Allele2, Allele, CAID) %>%
  filter(Drug %in% rownames(summary_df))
raw_anno <- raw_anno %>% rename(Category = PhenotypeCategory) %>% filter(Category != 'Other')
raw_anno[raw_anno$Category == 'Metabolism/PK',]$Category <- 'Metabolism'
length(unique(raw_anno$Gene))
length(unique(raw_anno$Drug))
length(unique(raw_anno$Variant))
length(unique(raw_anno$Allele))
length(unique(raw_anno$CAID))
# Colors
cat_color <- pal_npg("nrc")(9)[c(1,4,5,6)]
clinical_anno_df <- merge(clinical_anno_df, meta, on = 'Sample')
cat_order <- c('Toxicity', 'Dosage', 'Efficacy', 'Metabolism')#, 'Other')
raw_anno$Category <- factor(raw_anno$Category, levels = cat_order)

########### DRUG ###########
anno_count1 <- raw_anno[,c('Gene', 'Drug', 'Category')] %>% unique() %>%
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Gene')
anno_count2 <- raw_anno[,c('Gene', 'Variant', 'Drug', 'Category')] %>% unique() %>% #filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Variant')
anno_count3 <- raw_anno[,c('Gene', 'Allele', 'Drug', 'Category')] %>% unique() %>%#filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Allele')
anno_count <- rbind(anno_count1, anno_count2, anno_count3)
anno_count$Class <- factor(anno_count$Class, levels=c('Gene', 'Variant', 'Allele'))
#anno_count %>% ggplot(aes(x=Category, y=Count)) + geom_boxplot(aes(fill=Class)) +
#   scale_fill_manual('', values=cat_color) + theme_few() + 
#   labs(x='Phenotype category', y='Number of associations per drug')
anno_count %>% select(Category, Class, Count) %>% group_by(Category, Class) %>% 
  mutate(Mean = mean(Count)) %>% select(Category, Class, Mean) %>% unique()

########### Drug-x pairs ########### 
anno_count4 <- raw_anno[,c('Gene', 'Drug', 'Category', 'CAID')] %>% unique() %>%
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Gene x Annotation')
anno_count5 <- raw_anno[,c('Gene', 'Variant', 'Drug', 'Category', 'CAID')] %>% unique() %>% #filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Variant x Annotation')
anno_count6 <- raw_anno[,c('Gene', 'Allele', 'Drug', 'Category', 'CAID')] %>% unique() %>%#filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Allele x Annotation')

anno_count <- rbind(anno_count4, anno_count5, anno_count6)
anno_count$Class <- factor(anno_count$Class, levels=c('Gene x Annotation', 'Variant x Annotation', 'Allele x Annotation'))
anno_count %>% select(Category, Class, Count) %>% group_by(Category, Class) %>% 
  mutate(Mean = mean(Count)) %>% select(Category, Class, Mean) %>% unique()
# anno_count %>% ggplot(aes(x=Category, y=Count)) + geom_boxplot(aes(fill=Class)) +
#   scale_fill_manual('', values=cat_color) + theme_few() + 
#   labs(x='Phenotype categories', y='Number of genes or variants per drug')


########## Combine the above results
anno_count7 <- raw_anno[,c('Drug', 'Category', 'CAID')] %>% unique() %>%#filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Annotation')
library("scales")
show_col(pal_npg("nrc")(10))
ebtop <- function(x){return(mean(x)+plotrix::std.error(x)/sqrt(length(x)))}
ebbottom <- function(x){return(mean(x)-plotrix::std.error(x)/sqrt(length(x)))}

anno_count1 <- raw_anno[,c('Gene', 'Drug', 'Category')] %>% unique() %>%
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Gene')
anno_count2 <- raw_anno[,c('Variant', 'Drug', 'Category')] %>% unique() %>% #filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Variant')
anno_count4 <- raw_anno[,c('Gene', 'Drug', 'Category', 'CAID')] %>% unique() %>%
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Gene x Annotation')
anno_count5 <- raw_anno[,c('Variant', 'Drug', 'Category', 'CAID')] %>% unique() %>% #filter(!Gene %in% c('RYR1', 'CFTR', 'DPYD')) %>% 
  group_by(Category, Drug) %>% summarise(Count=n(), Class = 'Variant x Annotation')
anno_count <- rbind(anno_count1, anno_count2, anno_count4, anno_count5)#, anno_count7)
anno_count$Class <- factor(anno_count$Class, levels=c('Gene', 'Gene x Annotation', 'Variant', 'Variant x Annotation'))#, 'Annotation'))
anno_count %>% select(Category, Class, Count) %>% group_by(Category, Class) %>% 
  mutate(Mean = mean(Count)) %>% select(Category, Class, Mean) %>% unique()
pdf('./figure/Relation-G-V.pdf', width = 7, height = 3*0.85)
ggplot(anno_count, aes(x=Category, y=Count, fill=Class))+
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.8))+
  stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.8), width=0.3) +
  #scale_fill_brewer(palette = 'Paired') + theme_few() +
  scale_fill_manual('', values=c('#F39B7FFF', '#E64B35FF', '#8491B4FF', '#3C5488FF')) + #'#91D1C2FF', '#00A087FF')) +
  theme(legend.position = 'right') + theme_few() +
  labs(x='Phenotype category', y='Association count per drug')
dev.off()

# # Number of genes per drug
# p1 <- anno_count1 %>%
#   ggplot(aes(x=Category, y=count)) + 
#   geom_violin(aes(fill=Category))  + 
#   scale_fill_manual('', values=cat_color) + theme_few() + 
#   labs(x='Phenotype categories', y='Number of genes per drug')
# p2 <- anno_count2 %>%
#   ggplot(aes(x=Category, y=count)) + 
#   geom_violin(aes(fill=Category))  + 
#   scale_fill_manual('', values=cat_color) + theme_few() + 
#   labs(x='Phenotype categories', y='Number of variants per drug')
# p3 <- raw_anno[,c('Allele', 'Drug', 'Category')] %>% unique() %>%
#   group_by(Category, Drug) %>% summarise(count=n()) %>%
#   ggplot(aes(x=Category, y=count)) + 
#   geom_violin(aes(fill=Category))  + 
#   scale_fill_manual('', values=cat_color) + theme_few() + 
#   labs(x='Phenotype categories', y='Number of alleles per drug')


# ## CPIC 100 x 88
# library(ComplexHeatmap)
# colors <- c('Avoid'='#AF2318', 'Caution'='#F3B13E', 'Routine'='#2D695B', 'Missing'='grey')
# Heatmap(as.matrix(summary_df), na_col = "white", col = colors)

# res <- summary_df %>% mutate_all(funs(str_replace(., "Avoid", '-1'))) %>% 
#   mutate_all(funs(str_replace(., "Caution", '0'))) %>% 
#   mutate_all(funs(str_replace(., "Routine", '1'))) %>% 
#   mutate_all(funs(str_replace(., "Missing", '2')))
# colors <- c('-1'='#AF2318', '0'='#F3B13E', '1'='#2D695B', '2'='grey')
# Heatmap(as.matrix(res), na_col = "white", col = colors)


Avoid <- apply(summary_df,1,function(x) sum(x=='Avoid')/88)
Caution <- apply(summary_df,1,function(x) sum(x=='Caution')/88)
Routine <- apply(summary_df,1,function(x) sum(x=='Routine')/88)
Missing <- apply(summary_df,1,function(x) sum(x=='Missing')/88)
statis_df <- t(rbind(Routine, Caution, Avoid))
Heatmap(as.matrix(statis_df), na_col = "white", 
        #cluster_rows = FALSE, cluster_columns = FALSE,
        col = circlize::colorRamp2(seq(0, 1, length.out = 9), 
                                   colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(9)))






clinical_anno_df$PhenotypeCategory <- factor(clinical_anno_df$PhenotypeCategory, levels = cat_order)

anno_df <- clinical_anno_df %>% filter(PhenotypeCategory != 'Other')
anno_df_ratio <- anno_df %>% group_by(Sample) %>% 
  summarise(anno_count=n(), drug_count=length(unique(Drug))) %>% 
  mutate(ratio = anno_count/drug_count)
summary(anno_df_ratio$ratio)

# Number of clinical annotations for different samples
p1 <- anno_df %>% group_by(Group, Sample, PhenotypeCategory) %>% summarise(Count = n()) %>%
  ggplot(aes(x=Group, y=Count)) + geom_boxplot(aes(fill=Group), outlier.size=1) + 
  scale_fill_manual('', values=cat_color) + theme_few() + 
  labs(x='Biogeographical group', y='Number of clinical annotations') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_wrap(.~PhenotypeCategory, scales = 'free', ncol=5) + stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)], label.x.npc = 0.5)

# Number of drugs in different samples
p2 <- anno_df %>% select(Group, Sample, PhenotypeCategory, Drug) %>% distinct() %>%
  group_by(Group, Sample, PhenotypeCategory) %>% summarise(`Count` = n()) %>%
  ggplot(aes(x=Group, y=Count)) + geom_boxplot(aes(fill=Group), outlier.size=1)  + 
  scale_fill_manual('', values=cat_color) + theme_few() + 
  labs(x='Biogeographical group', y='Number of drugs') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_wrap(.~PhenotypeCategory, scales = 'free', ncol=5) + stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)], label.x.npc = 0.5)

pdf('figure5.pdf', width = 8, height = 8)
p1; p2; p2+theme(legend.position = 'right')
p + stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)], label.x.npc = 0.5)
dev.off()

##### Merge the above data
df1 <- anno_df %>% select(Group, Sample, PhenotypeCategory, Drug) %>% distinct() %>%
  group_by(Group, Sample, PhenotypeCategory) %>% summarise(`Drug` = n())
df2 <- anno_df %>% group_by(Group, Sample, PhenotypeCategory) %>% summarise(`Clinical Annotation` = n())
df <- merge(df1, df2)
df$Ratio <- df$`Clinical Annotation` / df$Drug
summary(df$Ratio)

colnames(df1) <- c('Group', 'Sample', 'PhenotypeCategory', 'Count'); colnames(df2) <- c('Group', 'Sample', 'PhenotypeCategory', 'Count')
df1$Class <- 'Drug'; df2$Class <- 'Clinical Annotation'
df_long <- rbind(df1, df2)#, df[c('Group', 'Sample', 'PhenotypeCategory', 'Ratio')])
p <- df_long %>% 
  ggplot(aes(x=Group, y=Count)) + geom_boxplot(aes(fill=Group), outlier.size=1)  + 
  scale_fill_manual('', values=cat_color) + theme_few() + theme_few() + 
  labs(x='Biogeographic group', y='Number per sample') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_grid(Class~PhenotypeCategory, scales = 'free_y')
p + stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)], label.x.npc = 0.5)

pdf('figure5a_merge.pdf', width = 8, height = 4)
p+theme(legend.position = 'bottom', legend.margin=margin(t = -0.2, unit='cm'))
dev.off()


#################################################################
# Number of genotypes in different samples
anno_df %>% group_by(Group, PhenotypeCategory, Sample, Alleles) %>% summarise(Count = n()) %>%
  ggplot(aes(x=Group, y=Count)) + geom_boxplot(aes(fill=Group)) + 
  scale_fill_manual('', values=cat_color) + theme_few() + 
  labs(x='Biogeographical group', y='Number of genotypes') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_wrap(.~PhenotypeCategory, scales = 'free', ncol=5) + 
  stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)])

# Number of "Genotype-Drug" in different samples: ignore this part.
anno_df['Genotype-Drug'] <- paste0(anno_df$Variant, anno_df$Alleles, anno_df$Drug)
anno_df %>% group_by(Group, PhenotypeCategory, Sample, `Genotype-Drug`) %>% distinct %>% summarise(Count = n()) %>%
  ggplot(aes(x=Group, y=Count)) + geom_boxplot(aes(fill=Group)) + 
  scale_fill_manual('', values=cat_color) + theme_few() + 
  labs(x='Biogeographical group', y='Number of drugs') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_wrap(.~PhenotypeCategory, scales = 'free', ncol=5)+#, space = 'free') + 
  stat_compare_means(aes(group=Group), label = "p.signif", fontface='bold', color=pal_npg("nrc")(9)[c(3)])


anno_df %>% ggdensity(x="EvidenceLevel", rug = TRUE, color="Group", size=1) + scale_color_viridis_d()



#################################################################
pgx_summary_df <- pgx_summary_df %>% filter(PhenotypeCategory != 'Other')
response_order <- c('Decreased', 'Moderate', 'Increased')
cat_order <- c('Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK')
pgx_summary_df$Response <- factor(pgx_summary_df$Response, levels = response_order)


## Number of occurrences of each drug in each category
drug_count <- pgx_summary_df %>% group_by(PhenotypeCategory, Drug) %>% summarise(Count = n())
# Counting by drug
drug_count %>% group_by(Drug) %>% summarise(Category.Num = n()) %>% filter(Category.Num >=3)

drug_freq <- pgx_summary_df %>% group_by(PhenotypeCategory, Drug) %>% summarise(Freq = n()/88)
drug_freq_ht <- drug_freq %>% pivot_wider(names_from = Drug, values_from = Freq) %>% as.data.frame()

row.names(drug_freq_ht) <- drug_freq_ht$PhenotypeCategory
drug_freq_ht <- drug_freq_ht[cat_order,][,-1]

pdf('figure5C.pdf', width = 20, height = 3.75)
Heatmap(as.matrix(drug_freq_ht), na_col = "white", 
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = circlize::colorRamp2(seq(0, 1, length.out = 9), colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(9)))
dev.off()
pdf('figure5C_slim.pdf', width = 8.5, height = 1.75)
input <- as.matrix(drug_freq_ht)
colnames(input) <- NULL
Heatmap(input, na_col = "white", 
        name = "Frequency",
        row_names_side = "left",
        column_title = paste(length(unique((drug_count$Drug))), 'drugs with high clinical annotation levels of evidence'), column_title_side = 'top',
        cluster_rows = FALSE, cluster_columns = FALSE,
        #col = circlize::colorRamp2(seq(0, 1, length.out = 11), colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))(11)))
        col = circlize::colorRamp2(seq(0, 1, length.out = 9), colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(9)))
dev.off()


#################################################################
#pgx_summary_df %>% group_by(PhenotypeCategory, EvidenceLevel, Response) %>% summarise(Freq = n())
#pgx_summary_df %>% group_by(PhenotypeCategory, Drug, EvidenceLevel) %>% summarise(Count = n())
count_detail <- pgx_summary_df %>% group_by(PhenotypeCategory, Drug, EvidenceLevel, Response) %>% summarise(count = n())

# Combine with drug_count to find the frequency of each drug in different categories
# in the cases where it already occurs
count_detail <- merge(count_detail, drug_count)
count_detail$Freq <- count_detail$count/count_detail$Count*100
# warfarin
count_detail[count_detail$Drug == 'warfarin',]

library(wesanderson)
names(wes_palettes)
purple <- c('#B6A0C3', '#6A4D94', '#413188')
count_detail$EvidenceLevel <- paste0('Level ', count_detail$EvidenceLevel)
cat_order <- c('Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK')
count_detail$PhenotypeCategory <- factor(count_detail$PhenotypeCategory, levels = cat_order)

g <- count_detail %>% 
  ggplot(aes(x=Response, y=Freq)) + geom_boxplot(aes(fill=Response), alpha=0.85, outlier.size=1) + 
  #scale_fill_brewer(palette = 'Paired') + scale_fill_manual(values = wes_palette("Royal2", n = 5)[3:5]) + 
  scale_fill_manual('', values = purple) + 
  theme_few() + 
  labs(x='Response Level', y='% Percentage') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom') + 
  facet_grid(EvidenceLevel~PhenotypeCategory, scales = 'free')

pdf('figure5e_merge.pdf', width = 8, height = 3.75)
g+theme(legend.position = 'bottom', legend.margin=margin(t = -0.2, unit='cm'))
dev.off()