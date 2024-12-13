#Introduction####
#Title: Founder Frequencies Visualization
#Purpose: Create supplementary figures that visualize the founder frequencies
#at QTL and not just the change in frequency. This script will create 
#supplemental QTL Frequencies figures. 
#Created: 4/10/23
#Last Edited: 6/23/23
#Packages Needed: 
library(tidyverse)
library(data.table)
library(ggpubr)
library(gridExtra)

#You do not need to create any additional directories as long as 
#"allhaps.zinc.015c.txt.gz" in your working directory. You will also need
#"ANOVA_haplotype_freq_scan_cM_QTLhits.txt" which can be made using Zinc Analysis
#R script. 

#Setup Founder Frequency Table#### 
input_haplo_freq_name <- "allhaps.zinc.015c.txt.gz"
haplo_freq = read.table(input_haplo_freq_name,header=TRUE)

haplo_freq_to_plot = haplo_freq %>%
  # GENERATE a TRT (treatment) vector that is simply a single
  # letter pulled from the replicate name
  # REMOVE 'pool' column
  mutate(TRT = substr(pool,1,1)) %>%
  select(-pool) %>%
  #     chr     pos    cM founder   freq TRT
  #1   chrX 1602150 0.924      A1 0.0267   C
  #2   chrX 1602150 0.924      A2 0.0031   C  
  
  # GROUP
  # CALCULATE MEAN of haplotype frequencies (over replicates)
  #   for each founder and marker interval
  dplyr::group_by(chr,cM,founder,TRT) %>%  
  summarize(freq = mean(freq)) #%>% Stop here for the frequency seperate 
# chr      cM founder TRT     freq
#   <chr> <dbl> <chr>   <chr>  <dbl>
#1 chr2L  2.61 A1      C     0.0216
#2 chr2L  2.61 A1      Z     0.0208

#Now have a table that has the founder's haplotype averaged frequency across all
#12 replicates

#Additional objects needed 
Founder_Color <- c("#F0E442", "#555555", "#E69F00", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

input_qtl_name <- "ANOVA_haplotype_freq_scan_cM_QTLhits.txt"
QTL = read.table(input_qtl_name,header=TRUE)
QTL$Label <- c("B", "C", "D", "E", "F", "G", "A")

#Frequency Side by Side Plot####
#Creates middle plot in the supplements
SbyS_plot_list = list()

for(i in 1:nrow(QTL)) {
  
  # ISOLATE just the QTL region in 'haplo_freq_to_plot'
  qtl_region = haplo_freq_to_plot %>% filter(chr==QTL$CHROM[i] &
                                                    cM > QTL$cM[i]-5 & cM < QTL$cM[i]+5)
  #Change C and Z to more informative names
  qtl_region$TRT[qtl_region$TRT == "C"] <- c("Control Non-Selected Pool")
  qtl_region$TRT[qtl_region$TRT == "Z"] <- c("Zinc Selected Pool")
  # GENERATE plot data
  
  SbyS_plot_list[[i]] = ggplot(qtl_region, aes(x=cM, y=freq, group=founder)) +
    facet_wrap(~TRT)+
    geom_smooth(data=qtl_region, method = "loess",
                aes(color=founder), se=FALSE, show.legend=TRUE) +
    scale_color_manual(values = Founder_Color)+
    ylab("Frequency") +
    xlab(paste0("Location (cM) on ",QTL$CHROM[i])) +
    # DRAW vertical black line showing where the QTL peak is
    geom_vline(xintercept = QTL$cM[i], linetype =
                 "dashed", colour = "black") +  
    # DRAW horizontal black line to mark where 0.125 is. (1/8 = 0.125)
    geom_hline(yintercept = 0.125, linetype = "dashed", colour="black")+
    labs(color="Founder")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_blank()) + 
    theme(panel.border = element_rect(fill = NA, color = "black"),
          axis.text=element_text(size=8),axis.title=element_text(size=10)) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent"))
}

#Run this if you want to see if all side by side plots
SbyS_plot_list

#If you want to see a particular QTL frequency
SbyS_plot_list[1]
#Change number 1-7 for each of QTL
#As a note: 1 corresponds to QTL B, and 7 corresponds to QTL A. 

#Frequency Peak Dot Plot####
#Creates the bottom plot of the supplements
freq_peak_list = list()

# LOOP through QTL
for(i in 1:nrow(QTL)) {
  
  # ISOLATE just the QTL region from 'haplo_freq_diff_to_plot'
  qtl_region = haplo_freq_to_plot %>% filter(chr==QTL$CHROM[i] &
                                                    cM == QTL$cM[i])
  qtl_region$TRT[qtl_region$TRT == "C"] <- c("Control Media")
  qtl_region$TRT[qtl_region$TRT == "Z"] <- c("25mM ZnCl2")
  # GENERATE plot data
  freq_peak_list[[i]] = ggplot(qtl_region, aes(x=founder, y=freq, group=TRT)) +
    # USE se=FALSE for lines without CIs on smoothing
    # USE se=TRUE  for lines with    CIs on smoothing (BUT
    #     will also get gray background on legend)
    #facet_wrap(~TRT)+
    geom_point(aes(color=TRT), size=2, position = position_dodge(0.1))+
    ylab("Frequency") +
    xlab("Founders") +
    ggtitle("Founder Frequencies at QTL Peak")+
    # DRAW vertical black line showing where the QTL peak is
    #geom_vline(xintercept = QTL$cM[i], linetype =
    #             "dashed", colour = "black") +  
    labs(color="Treatment")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_blank()) + 
    theme(panel.border = element_rect(fill = NA, color = "black"),
          axis.text=element_text(size=8),axis.title=element_text(size=10)) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent"))
}

#Run this if you want to see if all side by side plots
freq_peak_list

#If you want to see a particular QTL frequency
freq_peak_list[1]
#Change number 1-7 for each of QTL
#As a note: 1 corresponds to QTL B, and 7 corresponds to QTL A.

#Frequency Change Plot####
#Top plot of supplement.

# GET average haplotype frequencies (per founder per position)
# for each treatment (Selection / Control)
# THEN convert to a "Selection minus Control" difference
haplo_freq_diff_to_plot = haplo_freq %>%
  
  # GENERATE a TRT (treatment) vector that is simply a single
  #   letter pulled from the replicate name
  # REMOVE 'pool' column
  mutate(TRT = substr(pool,1,1)) %>%
  select(-pool) %>%
  
  # GROUP
  # CALCULATE MEAN of haplotype frequencies (over replicates)
  #   for each founder and marker interval
  group_by(chr,cM,founder,TRT) %>%
  summarize(freq = mean(freq)) %>%
  
  # REFORMAT to put the treatment-specific frequencies (per
  #   founder and position) in separate columns
  # CALCULATE DIFFERENCE (Selection minus Control)
  # REMOVE 'C' and 'Z' columns
  pivot_wider(names_from=TRT,values_from=freq) %>%
  mutate(Diff = Z-C) %>%
  select(-c(Z,C))

#End of pipeline we have ended up with differences in frequencies for each founder at each position 

Change_plot_list = list()

# LOOP through QTL
for(i in 1:nrow(QTL)) {
  
  # ISOLATE just the QTL region from 'haplo_freq_diff_to_plot'
  qtl_region = haplo_freq_diff_to_plot %>% filter(chr==QTL$CHROM[i] &
                                                    cM > QTL$cM[i]-5 & cM < QTL$cM[i]+5)
  #qtl_region$Diff <- qtl_region$Diff*100
  locpos = QTL$cM[i]
  
  # GENERATE plot data
  Change_plot_list[[i]] = ggplot(qtl_region, aes(x=cM, y=Diff, group=founder)) +
    # USE se=FALSE for lines without CIs on smoothing
    # USE se=TRUE  for lines with    CIs on smoothing (BUT
    #     will also get gray background on legend)
    scale_color_manual(values = Founder_Color)+
    geom_smooth(data=qtl_region, method = "loess",
               aes(color=founder), se=FALSE, show.legend=TRUE) +
    #geom_line(aes(color=founder), linewidth=1)+
    geom_vline(xintercept = locpos, linetype = "dashed", colour = "black") +
    ylab("Frequency Change") +
    xlab(paste0("Location (cM) on ",QTL$CHROM[i])) +
    labs(color="Founder")+
    theme_classic()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_blank()) + 
    theme(panel.border = element_rect(fill = NA, color = "black"),
          axis.text=element_text(size=8),axis.title=element_text(size=10)) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent"))
}

#Run this if you want to see if all side by side plots
Change_plot_list

#If you want to see a particular QTL frequency
Change_plot_list[1]
#Change number 1-7 for each of QTL
#As a note: 1 corresponds to QTL B, and 7 corresponds to QTL A.

#Making supplemental####
#So QTL A on X chromosome is plotted first. 
Numbers <- c(7, 1:6)

Graph_function_QTL <- function(CHANGE_PLOT, SBYS_PLOT, FREQ_PEAK, NUMBERS, QTL){
  pdf("Supplemental_2_QTL_Frequencies.pdf", width = 10, height = 10)
  for (i in NUMBERS) {
    #lay <- rbind(c(1,1),c(3,3),c(4,5))
    Graph_arrange <- grid.arrange(CHANGE_PLOT[[i]],
                                  SBYS_PLOT[[i]],
                                  FREQ_PEAK[[i]],
                                  ncol=1, top=text_grob(paste0("QTL ", QTL$Label[i])))
    Graph_arrange
  }
  graphics.off()
}

Graph_function_QTL(Change_plot_list, SbyS_plot_list, freq_peak_list, Numbers, QTL)
#Sometimes when running this function the page that should have QTL A is blank.
#If you rerun the function it usually will fix this. 
