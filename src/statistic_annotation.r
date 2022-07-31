library(ggplot2); library(ggsci); library(ggthemes); library(ggpubr); library(ComplexHeatmap); library(tidyverse)

clinical_anno_df <- read.csv('../res/clinical_anno_df.txt', sep="\t")
dosing_guideline_df <- read.csv('../res/dosing_guideline_df.txt', sep="\t")
pgx_summary_df <- read.csv('../res/pgx_summary_df.txt', sep="\t")

#race_dic = {'AFR': 'African American/Afro-Caribbean', 'AMR': 'Latino', 'EAS': 'East Asian', 'EUR': 'European'}
meta <- read.csv('../res/population.txt', sep='\t')
meta <- meta[c('Get.RM.137.Samples', 'Group')]; colnames(meta) <- c('Sample', 'Group')

cat_color <- pal_npg("nrc")(9)[c(1,4,5,6)]
clinical_anno_df <- merge(clinical_anno_df, meta, on = 'Sample')
cat_order <- c('Toxicity', 'Dosage', 'Efficacy', 'Metabolism/PK')#, 'Other')
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