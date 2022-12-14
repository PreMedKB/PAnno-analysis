library(ggplot2); library(ggsci); library(ggthemes); library(ggpubr); library(tidyverse)

######### Distribution
## Haplotype
haplotype <- data.frame()
race <- c("African American/Afro-Caribbean", "American", "Central/South Asian", "East Asian", "European", "Latino", "Near Eastern", "Oceanian", "Sub-Saharan African")


for (f in list.files('/Volumes/TM/PAnno/PAnno-db/data/pgx_diplotypes/population_frequency/')){
  gene <- strsplit(f, '_')[[1]][1]
  if (gene %in% c("CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
                  "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "NUDT15",
                  "SLCO1B1", "TPMT", "UGT1A1")){
    frequency <- read.csv(paste0('/Volumes/TM/PAnno/PAnno-db/data/pgx_diplotypes/population_frequency/', f), sep="\t")
    cols <- c("Allele", race)
    colnames(frequency) <- cols
    frequency$Gene <- gene
    haplotype <- rbind(haplotype, frequency)
  }
}


haplotype_frequency <- pivot_longer(haplotype, cols = -c(Allele, Gene), names_to = 'Group')
#haplotype_frequency$Allele <- factor(haplotype_frequency$Allele, levels = unique(haplotype_frequency$Allele))

## Diplotype
diplotype <- data.frame()
race <- c("African American/Afro-Caribbean", "American", "Central/South Asian", "East Asian", "European", "Latino", "Near Eastern", "Oceanian", "Sub-Saharan African")
for (f in list.files('/Volumes/TM/PAnno/PAnno-db/data/pgx_diplotypes/diplotype_frequency/')){
  gene <- strsplit(f, '_')[[1]][1]
  if (gene %in% c("CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
                  "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "NUDT15",
                  "SLCO1B1", "TPMT", "UGT1A1")){
    frequency <- read.csv(paste0('/Volumes/TM/PAnno/PAnno-db/data/pgx_diplotypes/diplotype_frequency/', f), sep="\t")
    cols <- c("Allele", race, "Global")
    colnames(frequency) <- cols
    frequency$Gene <- gene
    diplotype <- rbind(diplotype, frequency)
  }
}

diplotype_frequency <- pivot_longer(diplotype, cols = -c(Allele, Gene), names_to = 'Group')
#diplotype_frequency$Allele <- factor(diplotype_frequency$Allele, levels = unique(diplotype_frequency$Allele))


########################################
sd_per_allele <- diplotype_frequency %>% filter(Group != 'Global') %>% #value>2e-12, 
  group_by(Gene, Allele) %>% summarise(SD = sd(value), Mean = mean(value)) 
sd_per_allele$CV <- sd_per_allele$SD/sd_per_allele$Mean

# Boxplot of SD
tmp <- sd_per_allele %>% group_by(Gene) %>% summarise(median = median(SD)) %>% arrange(median)
sd_per_allele$Gene <- factor(sd_per_allele$Gene, levels = tmp$Gene)
p1 <- sd_per_allele %>%
  ggplot(aes(x=Gene, y=SD, color=Gene, fill=Gene, alpha=0.95)) +
  geom_boxplot() + theme_few() + 
  #theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) + 
  scale_color_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1) + 
  labs(y = 'SDs of diplotype frequencies', x=NULL)
p1+coord_flip()+theme(legend.position = "none")

# Boxplot of CV
#median_by_gene <- sd_per_allele %>% group_by(Gene) %>% summarise(median = median(CV)) %>% arrange(median)
#same_pattern <- read.csv('./same_pattern.txt', sep='\t')
#merge_df <- merge(median_by_gene, same_pattern)
#cor.test(merge_df$median, merge_df$Diplotype)


################### Figure 5a ######################
tmp <- sd_per_allele %>% group_by(Gene) %>% summarise(mean = mean(CV, na.rm = T)) %>% arrange(mean)
print(tmp)
sd_per_allele$Gene <- factor(sd_per_allele$Gene, levels = rev(tmp$Gene))
p2 <- sd_per_allele %>%
  ggplot(aes(x=Gene, y=CV, fill=Gene)) +
  geom_boxplot(alpha=0.9, outlier.size = 0.8) + theme_few() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) + 
  #scale_color_viridis_d(direction = 1) + 
  scale_fill_viridis_d(direction = 1) + 
  labs(y = 'CV of diplotype frequency', x=NULL)

pdf('/Volumes/TM/PAnno/PAnno-analysis/figure/CV_frequency_diplotype.pdf', width = 6.5, height = 2.7)
p2+theme(legend.position = "none")
dev.off()
#######################################################

# na.omit(sd_per_allele) %>% ggdensity(x="SD", rug = TRUE, color="Gene", size=1) + scale_color_viridis_d()## 这样看差异很小

cv1 <- sd_per_allele %>% group_by(Gene) %>% filter(CV>=1) %>% summarise(">=1" = n())
cv2 <- sd_per_allele %>% group_by(Gene) %>% filter(CV<1) %>% summarise("<1" = n()) 
cv_count_raw <- merge(cv1, cv2, all = TRUE)
cv_count <- pivot_longer(cv_count_raw, cols = -Gene, names_to = 'CV')
cv_count$CV <- factor(cv_count$CV, levels = c('>=1', '<1'))
p2 <- ggplot(cv_count, aes(x=Gene, y=value, fill=CV)) +
  geom_bar(position="stack", stat="identity", alpha=0.95) + 
  theme_few() + #scale_shape_manual("Group", values = shapes.race) + 
  #theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) + 
  scale_fill_nejm() + labs(y = 'Number of diplotypes', x=NULL) + scale_y_continuous(expand=c(0,0))
p2+coord_flip()


cv_count_raw[is.na(cv_count_raw)] <- 0
cv_count_raw$`CV>=1` <- cv_count_raw$`>=1`/(cv_count_raw$`>=1` + cv_count_raw$`<1`)*100
cv_count_raw <- cv_count_raw %>% arrange(`CV>=1`)
cv_count_raw$Gene <- factor(cv_count_raw$Gene, levels = cv_count_raw$Gene)
p4 <- ggplot(cv_count_raw, aes(x=Gene, y=`CV>=1`, fill=Gene)) +
  geom_bar(stat="identity", alpha=0.95) + 
  theme_few() + #theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) + 
  scale_fill_viridis_d(option = "plasma", direction = -1) + labs(y = 'Coefficient of Variation > 100%', x=NULL) + scale_y_continuous(expand=c(0,0))
p4

# Percentage
cv_count$Gene <- factor(cv_count$Gene, levels = cv_count_raw$Gene)
cv_count$CV <- factor(cv_count$CV, levels = c('<1', '>=1'))
p3 <- ggplot(cv_count, aes(x=Gene, y=value, fill=CV)) +
  geom_bar(position="fill", stat="identity", alpha=0.95, width=0.8) + 
  theme_few() + 
  #theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) + 
  scale_fill_brewer(palette = "Paired", direction = 1) + labs(y = 'Percentage of different CVs', x=NULL) + scale_y_continuous(expand=c(0,0))
p3#+coord_flip()


pdf('figure4_vertical_new.pdf', width = 4.5, height = 4.5)
p1+coord_flip()+theme(legend.position = "none")+theme(plot.margin=unit(c(1,5,1,1),'lines'))
p3+coord_flip()+theme(legend.position = "right")+theme(plot.margin=unit(c(1,1,1,1),'lines'))
dev.off()

pdf('figure4_long.pdf', width = 8, height = 3)
p1; p3+theme(legend.position = "none")
dev.off()

####################
shapes.race <- c(16, 18, 15, 17, 19, 8, 10, 13, 6); names(shapes.race) <- race
# Frequencies for haplotype and diplotye
p4 <- ggscatter(haplotype_frequency, "Allele", "value", 
                shape = 'Group',
                color = 'Group', alpha=0.9, palette = 'npg', 
                facet.by = 'Gene', scales = 'free', ncol = 4, 
                xlab = "Haplotype", ylab = "Frequency", legend.title = "Group") + 
  theme_few() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom") + 
  scale_shape_manual("Group", values = shapes.race)

p5 <- ggscatter(diplotype_frequency[diplotype_frequency$Group != 'Global',],
                "Allele", "value", 
                shape = 'Group',
                color = 'Group', alpha=0.9, palette = 'npg', 
                facet.by = 'Gene', scales = 'free', ncol = 4, 
                xlab = "Diplotype", ylab = "Frequency", legend.title = "Group") + 
  theme_few() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom") + 
  scale_shape_manual("Group", values = shapes.race)


pdf('supp_frequency.pdf', width = 10, height = 12.5)
p4; p5
dev.off()



# SD of each group
haplotype_frequency %>% group_by(Group, Gene) %>% summarise(Variance = sd(value)) %>%
  ggplot(aes(x=Gene, y=Variance)) +
  geom_point(aes(colour = Group, shape = Group), size = 3) + 
  theme_bw() + scale_shape_manual("Group", values = shapes.race) + 
  scale_color_npg() + labs(y = 'Standard deviation of haplotype frequencies')



#################### Figure 5b ####################
gene <- 'CYP4F2'
shapes.race <- c(16, 18, 15, 17, 19, 8, 10, 13, 6); names(shapes.race) <- race
h <- haplotype_frequency %>% filter(Gene == gene) %>% mutate(Class = 'Haplotypes')
d <- diplotype_frequency %>% filter(Gene == gene, Group != 'Global') %>% mutate(Class = 'Diplotypes')
hd <- rbind(h, d); hd$Class <- factor(hd$Class, levels = c('Haplotypes', 'Diplotypes'))
hd$`Biogeographic group` <-hd$Group
hdg <- ggplot(hd, aes(x=Allele, y=value)) +
  geom_point(aes(shape = `Biogeographic group`, color = `Biogeographic group`), 
             alpha = 1, size = 3, position = "jitter") + 
  theme_few() + scale_shape_manual(values = shapes.race) + 
  scale_color_npg() + labs(y = 'CYP4F2 frequency', x = NULL) + 
  facet_grid(.~Class, scales = 'free', space = 'free')


pdf('/Volumes/TM/PAnno/PAnno-analysis/figure/CYP4F2.pdf', width = 7.5*1.15, height = 3*0.9)
hdg #+ guides(color=guide_legend(nrow=3,bycol=TRUE)) + theme(legend.position = "bottom", legend.box="vertical") 
dev.off()
#######################################################

# indistinguishable diplotypes



