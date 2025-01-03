library(dplyr)  #1.1.4
library(ggplot2)  #3.4.4
library(tidyr)  #1.3.0
library(rstatix)  #0.7.2
library(ggpubr) #0.6.0
library(gridExtra)  #2.3
library(cowplot)  #1.1.1
library(ggforce)  #0.4.1
library(stringr)  #1.5.1

#import the data for analysis
all_reg <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/all_regions_merged_data", sep = " ", header = F)

colnames(all_reg) <- c("Region","Treat","Time","Sample","size","main_rep_size",
                       "nonscaled_cov","nonscaled_cutoff")


pOLA_cov <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/pOLA52_cov", sep = " ", header = F)
colnames(pOLA_cov) = c("Sample","pOLA52_cov")
chr_cov <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/chr_cov", sep = " ", header = F)
colnames(chr_cov) = c("Sample","chr_cov")


#Add annotations
inter_annot <- (read.csv("/Users/tkq300/OneDrive - University of Copenhagen/Users/tkq300/Desktop/mac_data/Mobilome_tool23/U00096_redone/all_regions_merged_data_annot_curated.csv", sep = ";"))
all_reg <- left_join(all_reg,inter_annot,by="Region")

#Remove IS1D. This is not the same as IS1F which is actually interesting. IS1D also occurs on pOLA52.
all_reg <- all_reg %>%
  filter(!coord == "1049828-1050546") #Coordinates of IS1D

#Fix annotations
all_reg$Annotation = gsub("^IS3 IS2 hyp. Prot$","IS3 IS2",
                          gsub("^IS3 IS3$", "IS3 IS3 ISEc17",
                               gsub("^tRNA$","tRNAs",
                                    gsub("NC_mhpE_mhpT","REP31a-g",
                                         gsub("NC_yaiL_frmB","REP32a-d",
                                              gsub("NC_tldD_yhdP","REP245a-g",
                                                   gsub("NC_rhaD_rhaA","REP299a-i",
                                                        gsub("NC_gltP_yjcO","REP321a-k",
                                                             gsub("NC_yjdN_yjdM","REP325a-i",
                                                                  gsub("NC_yjjV_yjjW","REP352a-h",
                                                                       gsub("NC_eco_mal","REP161a-i",
                                                                            gsub("ncRNA_ldrABC_rdlABC","ldr/rdl",
                                                                                 gsub("ncRNA_tRNA","tyrT operon",
                                                                                      gsub("tRNAs \\(serX\\)","tRNAs",
                                                                                           all_reg$Annotation))))))))))))))


all_reg$Annotation = gsub("^(IS1).$","\\1",
                          gsub("^(IS2).","\\1",
                               gsub("^(IS5).+","\\1",
                                    gsub("^(IS3).","\\1",
                                         gsub("IS10","IS150",
                                              all_reg$Annotation)))))


all_reg <- left_join(all_reg,pOLA_cov,by="Sample")
all_reg <- left_join(all_reg,chr_cov,by="Sample")

#Summarize and calculate exonuclease efficiency
all_reg.mean <- all_reg %>%
  group_by(Region,Sample,Time,Treat) %>%
  summarize(Mean_cov_increase = mean(as.numeric(nonscaled_cov)/as.numeric(nonscaled_cutoff)),
            Mean_cor_cov = mean(as.numeric(nonscaled_cov)-as.numeric(nonscaled_cutoff)),
            Annotation = Annotation,
            Relative_copies = (nonscaled_cov-nonscaled_cutoff)/pOLA52_cov,
            startcoord = gsub("-[0-9]*","",coord),
            Exonuclease_efficiency = (1-((chr_cov+1)/(pOLA52_cov/4+1)))*100
  )

all_reg.mean$Time <- factor(all_reg.mean$Time , levels = c("5","20","60"))

#Plot exonuclease efficiency. Supplementary figure 2.
all_reg.mean %>%
  group_by(Sample) %>%
  select(Sample,Exonuclease_efficiency) %>%
  unique() %>%
  ggplot(., aes(Sample, Exonuclease_efficiency)) +
  geom_col() +
  theme_minimal() +
  #ggtitle("Exonuclease efficiency") +
  ylab("Exonuclease efficiency") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylim(0,100)



all_reg.mean$Treat <- factor(all_reg.mean$Treat, levels = c("Cntr","Cet","Cu","Nal","SDS","Tet","UV"))

#Supplementary figure 3.
#All in one. With free y axis, separate e14 doesn't matter. 
all_reg.mean %>%
  #filter(!Annotation == "e14") %>%
  ggplot(.,aes(Treat,Mean_cov_increase,color = Time, type = Treat)) +
  geom_jitter(alpha = 1, width = 0.2, size = 2) + 
  facet_wrap(~Annotation, scales = "free_y") + 
  theme_bw() +
  geom_hline(yintercept = 1) +
  #scale_color_hue(direction = -1) + 
  scale_color_brewer(palette="Paired") +
  ylab("Mean increase in coverage\n over background") +
  xlab("Treatment")


all_treat_100int_roi <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/all_treat_100int_roi.bed", sep = "\t", header = F)
colnames(all_treat_100int_roi) = c("Treatment","Replicon","Start","End","Rep2","100_Start","100_End","Coverage","Region")

#Add annotations
all_treat_100int_roi.annot <- merge(all_treat_100int_roi,inter_annot, by="Region") 
all_treat_100int_roi.annot <- all_treat_100int_roi.annot %>%
  filter(!Annotation == "IS1D") #Remove IS1D


#Add other data
all_treat_100int_roi.annot$RegionTreat = paste(all_treat_100int_roi.annot$Region,all_treat_100int_roi.annot$Treatment, sep = "_")
all_reg$RegionTreat = paste(all_reg$Region,all_reg$Sample,sep = "_")
all_treat_100int_roi.annot.data <- merge(all_treat_100int_roi.annot, all_reg, by = "RegionTreat")
all_treat_100int_roi.annot.data <- merge(all_treat_100int_roi.annot.data, pOLA_cov, by.x = "Treatment", by.y = "Sample")

#Calculate relative ROI copies. Nonscaled cov and cutoff divided by pOLA coverage. Not adjusted for exo efficiency.
all_treat_100int_roi.annot.data$Relative_copies_ROI = (all_treat_100int_roi.annot.data$nonscaled_cov-all_treat_100int_roi.annot.data$nonscaled_cutoff)/all_treat_100int_roi.annot.data$pOLA52_cov.x

#Calculate exonuclease efficiency. It is a percentage efficiency. 100% means that chromosome is completely digested. 
#A pseudocount is added to both replicon coverages, to avoid division by 0.
all_treat_100int_roi.annot.data$Exo_efficiency = (1-((all_treat_100int_roi.annot.data$chr_cov+1)/(all_treat_100int_roi.annot.data$pOLA52_cov.x/4+1)))*100

#Replace negative exonuclease efficiency with 1 (they failed and should be given a low number)
temp <- all_treat_100int_roi.annot.data$Exo_efficiency
temp <- replace(temp, temp <1, 1)

all_treat_100int_roi.annot.data$Exo_efficiency = temp

#Get the total number of mapped bases
total_mapped_bases <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/total_mapped_bases", sep = " ", header = F)
colnames(total_mapped_bases) = c("Treatment","Total_mapped_bases")

#Add total mapped bases to data
all_treat_100int_roi.annot.data = left_join(all_treat_100int_roi.annot.data, total_mapped_bases, by = "Treatment")

#Adjust ROI increase in coverage by exonuclease efficiency. An efficiency of less than 100% means that increase over background (chr) is overestimated
all_treat_100int_roi.annot.data$Relative_copies_ROI_exo_efficiency = (all_treat_100int_roi.annot.data$Coverage-all_treat_100int_roi.annot.data$nonscaled_cutoff)*(all_treat_100int_roi.annot.data$Exo_efficiency/100)

#Correct pOLA52 coverage by exonuclease efficiency; then get ROI copy numbers per corrected pOLA52 coverage. 
#Both the pOLA52 and ROI coverage should be similarly affected by sequencing depth and both would be scaled by the same number
#Do not scale by pOLA52 copy number (4). This would leave out plasmids with unknown copy numbers and copy numbers varying according to growth stage. "Future studies should investigate these with e.g ddPCR".
all_treat_100int_roi.annot.data <- all_treat_100int_roi.annot.data %>%
  mutate(pOLA52_exo_corrected = pOLA52_cov.x * (Exo_efficiency/100), #Correct pOLA52 coverage for exo efficiency
         Rel_cop_pOLA52_exo_corrected = ((Relative_copies_ROI_exo_efficiency+1)/Total_mapped_bases)/
           ((pOLA52_exo_corrected+1)/Total_mapped_bases)) #Get ROI copies (adjusted for exo effic.) relative to adjusted pOLA52. Adjust both for total mapped bases (that cancels out)

#Remove Cet60, as it failed exonuclease treatment
all_treat_100int_roi.annot.data = all_treat_100int_roi.annot.data %>%
  filter(Treatment != "Cet60")

#Make table with basic stats
overall_table <- unique(all_treat_100int_roi.annot.data %>%
                          select(Treat,Time, Total_mapped_bases,chr_cov,pOLA52_cov.y,Exo_efficiency))
colnames(overall_table) <- c("Treatment","Time","Total mapped bases","Chr. cov.","pOLA52 cov.","Exonuclease eff.")

write.csv(overall_table,file = "/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/overall_table.csv")

png("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/overall_table.png")
overall_table_gt <- gt(overall_table, groupname_col = "Time", row_group_as_column = TRUE)
gtsave(overall_table_gt, "/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/overall_table.png")
dev.off

#Calculate means and replace <0 values (from subtracting the background cov)
Relative_copies_ROI_df <- all_treat_100int_roi.annot.data %>%
  select(Treatment,Annotation.y,Treat,Time,Coverage,nonscaled_cov,nonscaled_cutoff,
         Relative_copies_ROI,pOLA52_cov.x,chr_cov,Relative_copies_ROI_exo_efficiency,Rel_cop_pOLA52_exo_corrected,Exo_efficiency, coord.y) %>%
  group_by(Treatment, Annotation.y) %>%
  mutate(Relative_copies_ROI = mean(Relative_copies_ROI),
         Relative_copies_ROI_exo_efficiency = mean(Relative_copies_ROI_exo_efficiency),
         mean_Rel_cop_pOLA52_exo_corrected = mean(Rel_cop_pOLA52_exo_corrected)) %>%
  mutate(Relative_copies_ROI = replace(Relative_copies_ROI, which(Relative_copies_ROI<0), 0),
         Relative_copies_ROI_exo_efficiency = replace(Relative_copies_ROI_exo_efficiency, which(Relative_copies_ROI_exo_efficiency<0), 0),
         Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0),0)) #%>%

Relative_copies_ROI_df2 <- Relative_copies_ROI_df %>%
  select(Treatment,Annotation.y,Treat,Time,Relative_copies_ROI,
         pOLA52_cov.x,chr_cov,Coverage,nonscaled_cov,nonscaled_cutoff,coord.y) %>%
  group_by(Treatment, Annotation.y) %>%
  mutate(mean_rel_copies = mean(Relative_copies_ROI),
         simple_increase = Coverage/nonscaled_cutoff) %>%
  select(Treatment, Annotation.y,Treat,Time,mean_rel_copies,pOLA52_cov.x,chr_cov,Coverage,nonscaled_cov,nonscaled_cutoff,simple_increase,Relative_copies_ROI,coord.y) %>%
  distinct()

Relative_copies_ROI_df3 <- Relative_copies_ROI_df2 %>%
  select(Treatment,Annotation.y,Treat,Time,Relative_copies_ROI,mean_rel_copies,pOLA52_cov.x,chr_cov) %>%
  distinct()

Relative_copies_ROI_df3$Treat <- factor(Relative_copies_ROI_df3$Treat, 
                                        levels = c("Cntr","Cet","Cu","Nal","SDS","Tet","UV"))


Relative_copies_ROI_df4 <- Relative_copies_ROI_df %>%
  select(Treat, Time, Annotation.y, Rel_cop_pOLA52_exo_corrected,coord.y) %>%
  group_by(Treat,Time,Annotation.y,coord.y) %>%
  mutate(mean_Rel_cop_pOLA52_exo_corrected = mean(Rel_cop_pOLA52_exo_corrected)) %>%
  select(Treat, Time, Annotation.y, mean_Rel_cop_pOLA52_exo_corrected,coord.y) %>%
  distinct()
Relative_copies_ROI_df4$Treat <- factor(Relative_copies_ROI_df4$Treat, 
                                        levels = c("Cntr","Cet","Cu","Nal","SDS","Tet","UV"))

#Get means for all ROIs within groups
Relative_copies_ROI_df5 <- Relative_copies_ROI_df %>%
  select(Treat, Time, Annotation.y, Rel_cop_pOLA52_exo_corrected) %>%
  group_by(Treat,Time,Annotation.y) %>%
  mutate(mean_Rel_cop_pOLA52_exo_corrected = mean(Rel_cop_pOLA52_exo_corrected)) %>%
  select(Treat, Time, Annotation.y, mean_Rel_cop_pOLA52_exo_corrected) %>%
  distinct()

#Figure 1.
ggplot(Relative_copies_ROI_df4, aes(x=factor(Time, level = c("5","20","60")),mean_Rel_cop_pOLA52_exo_corrected, color = Treat)) +
  geom_jitter(alpha = 1, width = 0.2, size = 1.5) + 
  facet_wrap(~Annotation.y,scales = "free_y") +
  theme_bw() +
  xlab("Time") +
  scale_color_brewer(palette="Paired") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                     name = "ROI copies relative to pOLA52\nadjustedd for exonuclease efficiency")

all_treat_100int_roi.annot.data$Treat <- factor(all_treat_100int_roi.annot.data$Treat, 
                                                levels = c("Cntr","Cet","Cu","Nal","SDS","Tet","UV"))

all_treat_100int_roi.annot.data$Corrected_increase = 
  all_treat_100int_roi.annot.data$nonscaled_cov/all_treat_100int_roi.annot.data$nonscaled_cutoff


all_treat_100int_roi.annot.data$Relative_copies = (all_treat_100int_roi.annot.data$nonscaled_cov-all_treat_100int_roi.annot.data$nonscaled_cutoff)/all_treat_100int_roi.annot.data$pOLA52_cov.x

#Check the mean increase in scaled coverage
print(all_treat_100int_roi.annot.data %>%
        group_by(Annotation.y, Treatment) %>%
        summarise(mean_increase = sum(Corrected_increase)) %>%
        arrange(desc(mean_increase)), n = 100)

#e14 and Stress region_2 are the biggest changes (they are also largest)! Followed by rRNA. 
#Check regions for filtering
summed_changes <- all_treat_100int_roi.annot.data %>%
  group_by(Annotation.y) %>%
  summarise(median = median(Corrected_increase),
            sd = sd(Corrected_increase),
            mean = mean(Corrected_increase),
            median_sd = median(Corrected_increase) + sd(Corrected_increase),
            sum = sum(Corrected_increase),
            max = max(Corrected_increase),
            copies = mean(Relative_copies)) %>%
  arrange(desc(sum))
print(summed_changes, n=100)
#Stress region_1 and Hyp. Prot_1 are very low

#Filter out very high changing regions and then plot density. Try to figure out where the cutoff should be.

#summed_changes_filt <- summed_changes %>%
#  filter(median > 1 | sd > 1)
#Removes many that are either not varying much or have low median change

ggplot(summed_changes, aes(max)) +
  geom_density() +
  geom_vline(xintercept = 3)



#### NOT USED REALLY###
summed_changes_filt <- summed_changes %>%
  #filter(max > 2)
  filter(max > 0)


#Fix annotations
summed_changes_filt$Annotation.y = gsub("^IS3 IS2 hyp. Prot$","IS3 IS2",
                                        gsub("^IS3 IS3$", "IS3 IS3 ISEc17",
                                             gsub("^tRNA$","tRNAs",
                                                  gsub("NC_mhpE_mhpT","REP31a-g",
                                                       gsub("NC_yaiL_frmB","REP32a-d",
                                                            gsub("NC_tldD_yhdP","REP245a-g",
                                                                 gsub("NC_rhaD_rhaA","REP299a-i",
                                                                      gsub("NC_gltP_yjcO","REP321a-k",
                                                                           gsub("NC_yjdN_yjdM","REP325a-i",
                                                                                gsub("NC_yjjV_yjjW","REP352a-h",
                                                                                     gsub("NC_eco_mal","REP161a-i",
                                                                                          gsub("ncRNA_ldrABC_rdlABC","ldr/rdl",
                                                                                               gsub("ncRNA_tRNA","tyrT operon",
                                                                                                    gsub("tRNAs \\(serX\\)","tRNAs",
                                                                                                         summed_changes_filt$Annotation.y))))))))))))))


regions_keep = summed_changes_filt$Annotation.y


#The filtering was changed and currently no regions are removed based on filtering.
all_treat_100int_roi.annot.data.filt = all_treat_100int_roi.annot.data %>%
  filter(Annotation.y %in% regions_keep)

all_treat_100int_roi.annot.data.filt$Annotation.y <- factor(all_treat_100int_roi.annot.data.filt$Annotation.y , levels = regions_keep)

#Split up the regions
is_elements <- c("IS1","IS150","IS2","IS3","IS5") #Exclude Hyp. Prot_2 as it is the same as mcrB (which should be mrcB) 
misc_elements <- c("rhsA","rhsB","rhsC","rhsD","rhsE") 
REP_elements <- c("REP31a-g","REP325a-i","REP161a-i","REP32a-d","REP299a-i","REP321a-k","REP352a-h","REP245a-g")
e14_element <- "e14"
stress_regions <- c("ldr/rdl","tyrT operon","rRNA","tRNAs")

#Replace negative values with 0 and then add tiny pseudocount. 
all_treat_100int_roi.IS <- all_treat_100int_roi.annot.data.filt %>% 
  filter(Annotation.y %in% is_elements) %>%
  mutate(Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0), 0.0))

all_treat_100int_roi.misc <- all_treat_100int_roi.annot.data.filt %>%
  filter(Annotation.y %in% misc_elements) %>%
  mutate(Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0), 0.0))

all_treat_100int_roi.REP <- all_treat_100int_roi.annot.data.filt %>%
  filter(Annotation.y %in% REP_elements) %>%
  mutate(Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0), 0.0))

#Rename REP elements
all_treat_100int_roi.REP$Annotation.y = gsub("NC_mhpE_mhpT","REP31a-g",
                                             gsub("NC_yaiL_frmB","REP32a-d",
                                                  gsub("NC_eco_mal","REP161a-i",
                                                       gsub("NC_tldD_yhdP","REP245a-g",
                                                            gsub("NC_rhaD_rhaA","REP299a-i",
                                                                 gsub("NC_gltP_yjcO","REP321a-k",
                                                                      gsub("NC_yjdN_yjdM","REP325a-i",
                                                                           gsub("NC_yjjV_yjjW","REP352a-h",
                                                                                all_treat_100int_roi.REP$Annotation.y))))))))

all_treat_100int_roi.e14 <- all_treat_100int_roi.annot.data.filt %>%
  filter(Annotation.y %in% e14_element) %>%
  mutate(Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0), 0.0))

all_treat_100int_roi.stress <- all_treat_100int_roi.annot.data.filt %>%
  filter(Annotation.y %in% stress_regions) %>%
  mutate(Rel_cop_pOLA52_exo_corrected = replace(Rel_cop_pOLA52_exo_corrected, which(Rel_cop_pOLA52_exo_corrected<0), 0.0))

#850x650 export as svg
#Figure 3.
ggboxplot(all_treat_100int_roi.IS, x = "Treat", y = "Rel_cop_pOLA52_exo_corrected", color = "Treat", outlier.shape = NA) +
  rotate_x_text(angle = 45) +
  geom_jitter(shape = 16, alpha = 0.3, size = 0.5) +
  facet_grid(Annotation.y~Time,scales = "free_y") +
  stat_summary(fun=mean, geom="point", shape=17, size=1) +
  theme_bw() +
  scale_color_brewer(palette="Paired") +
  scale_y_sqrt(limits = c(-0.05,max(all_treat_100int_roi.IS$Rel_cop_pOLA52_exo_corrected)),name = "ROI/pOLA52 adjusted for exonuclease efficiency\n (square-root transformed)") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top",
                     label.x.npc = "left",
                     vjust = 1,
                     #hjust = -0.1,
                     label = "p.format") +
  stat_pwc(label = "p.adj.signif", method = "wilcox.test",
           #method.args = list(exact = "TRUE"),
           ref.group = "Cntr", 
           y.position = -0.05,
           #y.position = y_pos,
           #vjust = 1,
           hide.ns = T,
           p.adjust.method = "bonferroni", 
           remove.bracket = T,
           step.increase = 0,
           dodge = 0,
           group.by = "legend.var",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("a", "b", "c", "d", "ns")))

temp <- all_treat_100int_roi.IS %>% filter(Annotation.y == "IS1" & Time == "60")
pairwise.wilcox.test(sqrt(temp$Rel_cop_pOLA52_exo_corrected), temp$Sample, p.adjust.method="bonferroni",exact = F)

#Figure 2.
ggboxplot(all_treat_100int_roi.e14, x = "Treat", y = "Rel_cop_pOLA52_exo_corrected", color = "Treat", outlier.shape = NA) +
  rotate_x_text(angle = 45) +
  geom_jitter(shape = 16, alpha = 0.3, size = 0.5) +
  rotate_x_text(angle = 45) +
  facet_grid(Annotation.y~Time,scales = "free_y") +
  stat_summary(fun=mean, geom="point", shape=17, size=1) +
  theme_bw() +
  scale_color_brewer(palette="Paired") +
  scale_y_sqrt(limits = c(-0.02,max(all_treat_100int_roi.e14$Rel_cop_pOLA52_exo_corrected)),name = "ROI/pOLA52 adjusted for exonuclease efficiency\n (square-root transformed)") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top",
                     label.x.npc = "left",
                     vjust = 1,
                     #hjust = -0.1,
                     label = "p.format") +
  stat_pwc(label = "p.adj.signif", method = "wilcox.test",
           #method.args = list(exact = "TRUE"),
           ref.group = "Cntr", 
           y.position = -0.045,
           #y.position = y_pos,
           #vjust = 1,
           hide.ns = T,
           p.adjust.method = "bonferroni", 
           remove.bracket = T,
           step.increase = 0,
           dodge = 0,
           group.by = "legend.var",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("a", "b", "c", "d", "ns")))

temp <- all_treat_100int_roi.e14 %>% filter(Annotation.y == "e14" & Time == "20")
pairwise.wilcox.test(sqrt(temp$Rel_cop_pOLA52_exo_corrected), temp$Sample, p.adjust.method="bonferroni",exact = F)

#Figure 4.
ggboxplot(all_treat_100int_roi.REP, x = "Treat", y = "Rel_cop_pOLA52_exo_corrected", color = "Treat", outlier.shape = NA) +
  rotate_x_text(angle = 45) +
  geom_jitter(shape = 16, alpha = 0.3, size = 0.5) +
  rotate_x_text(angle = 45) +
  facet_grid(Annotation.y~Time,scales = "free_y") +
  stat_summary(fun=mean, geom="point", shape=17, size=1) +
  theme_bw() +
  scale_color_brewer(palette="Paired") +
  scale_y_sqrt(limits = c(-0.05,max(all_treat_100int_roi.REP$Rel_cop_pOLA52_exo_corrected)),name = "ROI/pOLA52 adjusted for exonuclease efficiency\n (square-root transformed)") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top",
                     label.x.npc = "left",
                     vjust = 1,
                     #hjust = -0.1,
                     label = "p.format") +
  stat_pwc(label = "p.adj.signif", method = "wilcox.test",
           #method.args = list(exact = "TRUE"),
           ref.group = "Cntr", 
           y.position = -0.05,
           #y.position = y_pos,
           #vjust = 1,
           hide.ns = T,
           p.adjust.method = "bonferroni", 
           remove.bracket = T,
           step.increase = 0,
           dodge = 0,
           group.by = "legend.var",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("a", "b", "c", "d", "ns")))


temp <- all_treat_100int_roi.REP %>% filter(Annotation.y == "REP245a-g" & Time == "60")
pairwise.wilcox.test(sqrt(temp$Rel_cop_pOLA52_exo_corrected), temp$Sample, p.adjust.method="bonferroni",exact = F)

#Supplementary figure 8.
ggboxplot(all_treat_100int_roi.stress, x = "Treat", y = "Rel_cop_pOLA52_exo_corrected", color = "Treat", outlier.shape = NA) +
  rotate_x_text(angle = 45) +
  geom_jitter(shape = 16, alpha = 0.3, size = 0.5) +
  rotate_x_text(angle = 45) +
  facet_grid(Annotation.y~Time,scales = "free_y") +
  stat_summary(fun=mean, geom="point", shape=17, size=1) +
  theme_bw() +
  scale_color_brewer(palette="Paired") +
  scale_y_sqrt(limits = c(-0.05,max(all_treat_100int_roi.stress$Rel_cop_pOLA52_exo_corrected)),name = "ROI/pOLA52 adjusted for exonuclease efficiency\n (square-root transformed)") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top",
                     label.x.npc = "left",
                     vjust = 1,
                     #hjust = -0.1,
                     label = "p.format") +
  stat_pwc(label = "p.adj.signif", method = "wilcox.test",
           #method.args = list(exact = "TRUE"),
           ref.group = "Cntr", 
           y.position = -0.05,
           #y.position = y_pos,
           #vjust = 1,
           hide.ns = T,
           p.adjust.method = "bonferroni", 
           remove.bracket = T,
           step.increase = 0,
           dodge = 0,
           group.by = "legend.var",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("a", "b", "c", "d", "ns")))

#Supplementary figure 11.
ggboxplot(all_treat_100int_roi.misc, x = "Treat", y = "Rel_cop_pOLA52_exo_corrected", color = "Treat", outlier.shape = NA) +
  rotate_x_text(angle = 45) +
  geom_jitter(shape = 16, alpha = 0.3, size = 0.5) +
  rotate_x_text(angle = 45) +
  facet_grid(Annotation.y~Time,scales = "free_y") +
  stat_summary(fun=mean, geom="point", shape=17, size=1) +
  theme_bw() +
  scale_color_brewer(palette="Paired") +
  scale_y_sqrt(limits = c(-0.05,max(all_treat_100int_roi.misc$Rel_cop_pOLA52_exo_corrected)),name = "ROI/pOLA52 adjusted for exonuclease efficiency\n (square-root transformed)") +
  stat_compare_means(method = "kruskal.test", 
                     label.y.npc = "top",
                     label.x.npc = "left",
                     vjust = 1,
                     #hjust = -0.1,
                     label = "p.format") +
  stat_pwc(label = "p.adj.signif", method = "wilcox.test",
           #method.args = list(exact = "TRUE"),
           ref.group = "Cntr", 
           y.position = -0.05,
           #y.position = y_pos,
           #vjust = 1,
           hide.ns = T,
           p.adjust.method = "bonferroni", 
           remove.bracket = T,
           step.increase = 0,
           dodge = 0,
           group.by = "legend.var",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("a", "b", "c", "d", "ns")))

#For a few ROIs, tests do not work. Here is the stats and sig levels are added manually to figures in Inkscape.
#Remember to change time to 5, 20, and 60 and get the individual stats.
#Kruskal-Wallis test is insignificant, no pairwise tests performed (e.g. REP299a-i)
temp <- all_treat_100int_roi.e14 %>% filter(Annotation.y == "e14" & Time == "60")
pairwise.wilcox.test(temp$Rel_cop_pOLA52_exo_corrected, temp$Sample, p.adjust.method="bonferroni",exact = F)


temp <- all_treat_100int_roi.stress %>% filter(Annotation.y == "tyrT operon" & Time == "5")
pairwise.wilcox.test(temp$Rel_cop_pOLA52_exo_corrected, temp$Sample, p.adjust.method="bonferroni",exact = F)

temp <- all_treat_100int_roi.misc %>% filter(Annotation.y == "rhsE" & Time == "60")
pairwise.wilcox.test(temp$Rel_cop_pOLA52_exo_corrected, temp$Sample, p.adjust.method="bonferroni",exact = F)



#e14 copies in UV60
all_treat_100int_roi.e14 %>%
  filter(Treatment == "UV60") %>%
  summarize((Rel_cop_pOLA52_exo_corrected))

####Investigate clipped and discordant reads within ROIs
###Should the discordant data be normalized somehow? E.g. per million total mapped bases?
#Clipped reads
all_clip <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/clip_cov_roi_all", sep = "\t", header = F) 
#Discordant reads
all_clip <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/disc_cov_roi_all", sep = "\t", header = F) 
colnames(all_clip) = c("Replicon","Position","Clipped_coverage","Treat","Time","Sample","Region","Terminus")

#Control regions
all_clip_cntr <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/disc_cov_roi_all_cntr", sep = "\t", header = F) 
colnames(all_clip_cntr) = c("Replicon","Position","Clipped_coverage","Treat","Time","Sample","Region","Terminus")


cntr_regions <- c("U00096.3:4181245-4185273","U00096.3:4056625-4058034","U00096.3:1765629-1768685",
                  "U00096.3:1632277-1652745","U00096.3:68348-70048")
cntr_annot <- c("rpoB","glnA","ydiJ","Qin","araB")
cntr_annot <- cbind(as.data.frame(cntr_regions),
                    as.data.frame(cntr_annot)
)
colnames(cntr_annot) <- c("Region","Annotation")


total_mapped_bases <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/total_mapped_bases", sep = " ", header = F)
total_mapped_clip_bases <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/total_mapped_clip_bases", sep = " ", header = F)
total_mapped_clip_bases <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/total_mapped_disc_bases", sep = " ", header = F)

colnames(total_mapped_clip_bases) = c("Sample","total_mapped_clip_bases") 
colnames(total_mapped_bases) = c("Sample","total_mapped_bases")

both_mapped_bases <- merge(total_mapped_clip_bases, total_mapped_bases, by = "Sample")
both_mapped_bases$relative_clipped <- both_mapped_bases$total_mapped_clip_bases/both_mapped_bases$total_mapped_bases*100

all_clip.annot <- left_join(all_clip,inter_annot,by="Region")
all_clip.annot <- left_join(all_clip.annot,both_mapped_bases,by="Sample")
all_clip.annot <- left_join(all_clip.annot,all_reg, by = c("Region","Sample"))


#Add nonscaled_cov for given regions and scale the clipped coverage with this. 

#Control regions processing
all_cntr <- read.csv("/Users/tkq300/mac_data/Mobilome_tool23/U00096_redone/all_regions_merged_data_cntr", sep = " ", header = F)

#For U00096, the replicon size was not included
#all_reg$main_rep_size = 4641652

colnames(all_cntr) <- c("Region","Treat","Time","Sample","size","main_rep_size",
                        "nonscaled_cov")
all_clip_cntr.annot <- left_join(all_clip_cntr, all_cntr, by=c("Region","Sample"))
all_clip_cntr.annot <- left_join(all_clip_cntr.annot,cntr_annot,by="Region")


#Fix annotations
all_clip.annot$Annotation.x = gsub("^(IS1).","\\1",
                                   gsub("^(IS2).","\\1",
                                        gsub("^(IS5).+","\\1",
                                             gsub("^(IS3).","\\1",
                                                  all_clip.annot$Annotation.x))))

all_clip.annot$Annotation.x = gsub("^IS3 IS2 hyp. Prot$","IS3 IS2",
                                   gsub("^IS3 IS3$", "IS3 IS3 ISEc17",
                                        gsub("^tRNA$","tRNAs",
                                             gsub("NC_mhpE_mhpT","REP31a-g",
                                                  gsub("NC_yaiL_frmB","REP32a-d",
                                                       gsub("NC_tldD_yhdP","REP245a-g",
                                                            gsub("NC_rhaD_rhaA","REP299a-i",
                                                                 gsub("NC_gltP_yjcO","REP321a-k",
                                                                      gsub("NC_yjdN_yjdM","REP325a-i",
                                                                           gsub("NC_yjjV_yjjW","REP352a-h",
                                                                                gsub("NC_eco_mal","REP161a-i",
                                                                                     gsub("ncRNA_ldrABC_rdlABC","ldr/rdl",
                                                                                          gsub("ncRNA_tRNA","tyrT operon",
                                                                                               gsub("tRNAs \\(serX\\)","tRNAs",
                                                                                                    all_clip.annot$Annotation.x))))))))))))))

#Get the 100 bases around the start and end position of each ROI
all_clip.annot = all_clip.annot %>%
  group_by(Region,Sample,Terminus) %>%
  mutate(Index = row_number()-100)

all_clip_cntr.annot = all_clip_cntr.annot %>%
  group_by(Region,Sample,Terminus) %>%
  mutate(Index = row_number()-100)

#Scale clipped coverage by the total coverage
all_clip.annot = all_clip.annot %>%
  group_by(Region,Sample,Terminus) %>%
  mutate(Corrected_clip = Clipped_coverage/nonscaled_cov)

all_clip_cntr.annot = all_clip_cntr.annot %>%
  group_by(Region,Sample,Terminus) %>%
  mutate(Corrected_clip = Clipped_coverage/nonscaled_cov)

unique(all_clip.annot$Annotation.x)

#Rename "start" and "end" to left and right for Terminus
all_clip.annot$Terminus = gsub("start","Left",
                               gsub("end","Right", all_clip.annot$Terminus))
is_elements <- c("IS1","IS150","IS2","IS3","IS5") #Exclude Hyp. Prot_2 as it is the same as mcrB (which should be mrcB) 
misc_elements <- c("rhsA","rhsB","rhsC","rhsD","rhsE") 
REP_elements <- c("REP31a-g","REP325a-i","REP161a-i","REP32a-d","REP299a-i","REP321a-k","REP352a-h","REP245a-g")
e14_element <- "e14"
stress_regions <- c("ldr/rdl","tyrT operon","rRNA","tRNAs")

#Replace Corrected_clip values > 1 with 1. These are due to VERY low average nonscaled_cov values.
all_clip.annot$Corrected_clip[all_clip.annot$Corrected_clip > 1] <- 1

all_clip.annot.IS <- all_clip.annot %>% 
  filter(Annotation.x %in% is_elements)

all_clip.annot.misc <- all_clip.annot %>%
  filter(Annotation.x %in% misc_elements)

all_clip.annot.REP <- all_clip.annot %>%
  filter(Annotation.x %in% REP_elements)

all_clip.annot.e14 <- all_clip.annot %>%
  filter(Annotation.x %in% e14_element)

all_clip.annot.stress <- all_clip.annot %>%
  filter(Annotation.x %in% stress_regions)

####Control regions - it is now disc and not clipped reads.
#Supplementary figure 1.
ggplot(all_clip_cntr.annot, aes(Index,Corrected_clip, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation~Sample) +
  scale_x_continuous(breaks = c(-100,0,100)) +
  theme_bw() +
  rotate_x_text(angle = 90) +
  ylab("Corrected discordant reads") +
  xlab("Index position")

#Supplementary figure 6.
ggplot(all_clip.annot.IS, aes(Index,Corrected_clip, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation.x~Sample,scales = "free_y") +
  scale_x_continuous(breaks = c(-100,0,100)) +
  theme_bw() +
  rotate_x_text(angle = 90) +
  ylab("Corrected discordant reads") +
  xlab("Index position")

#Supplementary figure 12.
ggplot(all_clip.annot.misc, aes(Index,Corrected_clip, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation.x~Sample,scales = "free_y") + 
  theme_bw() +
  scale_x_continuous(breaks = c(-100,0,100)) +
  rotate_x_text(angle = 90) +
  ylab("Corrected discordant reads") +
  xlab("Index position")

#Supplementary figure 7.
ggplot(all_clip.annot.REP, aes(Index,Corrected_clip, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation.x~Sample,scales = "free_y") + 
  theme_bw() +
  scale_x_continuous(breaks = c(-100,0,100)) +
  rotate_x_text(angle = 90) +
  ylab("Corrected discordant reads") +
  xlab("Index position")

#Supplementary figure 8.
ggplot(all_clip.annot.stress, aes(Index,Corrected_clip, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation.x~Sample,scales = "free_y") + 
  theme_bw() +
  scale_x_continuous(breaks = c(-100,0,100)) +
  rotate_x_text(angle = 90) +
  ylab("Corrected discordant reads") +
  xlab("Index position")

#Not included in manuscript
ggplot(all_clip.annot.e14, aes(Index,Clipped_coverage, color = Terminus)) +
  geom_smooth(method = "loess") +
  facet_grid(Annotation.x~Sample) +
  theme_bw() +
  scale_x_continuous(breaks = c(-100,0,100)) +
  rotate_x_text(angle = 90)