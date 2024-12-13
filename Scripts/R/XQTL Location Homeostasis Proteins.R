#Introduction####
#Title: Location of Homeostasis Proteins relative to QTL
#Purpose: Look where known zinc homeostasis genes (ZnTs, ZIPs, MTF-1, Mtn) are 
#in comparison to our QTL and the manually filtered our peaks. 
#Create: 6/5/23
#Last Edited: 6/23/23
#Packages Needed:
library(tidyverse)

#You will need the "Known Zinc homeostasis Genes.txt" 
#You will also need "ANOVA_haplotype_freq_scan_cM_QTLhits.txt" and 
#ANOVA_haplotype_freq_scan_cM_smoo.txt" which were made with Zinc Analysis
#R Code.
#Additionally, you will need the "ANOVA_haplotype_freq_scan_cM_Filtered_Hits.txt"
#which was the optional step in Zinc Analysis R Code at line 329

#Load in data sets####
homeo_proteins <- read.table("Known Zinc homeostasis Genes.txt", header = TRUE)
homeo_proteins$Chromosome_Recomb <- as.character(homeo_proteins$Chromosome_Recomb)
#This table and its information was manually collected from Flybase. 
# version FB2023_01, released February 15, 2023
head(homeo_proteins)
#     Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 1 ZnT33D   ZnT chr2L   12164241    12166704                 2     46
# 2 ZnT35C   ZnT chr2L   15228728    15247714                 2     51
# 3 ZnT41F   ZnT chr2R    5774776     5782858                 2     55
# 4 ZnT49B   ZnT chr2R   12652805    12655743                 2     67
# 5 ZnT63C   ZnT chr3L    3303708     3314009                 3      7
# 6 ZnT77C   ZnT chr3L   20484293    20494921                 3     47

Genomic_smoothed <- read.table("ANOVA_haplotype_freq_scan_cM_smoo.txt", header = TRUE)
QTL_cM <- read.table("ANOVA_haplotype_freq_scan_cM_QTLhits.txt", header = TRUE)
Filt_Peaks_cM <- read.table("ANOVA_haplotype_freq_scan_cM_Filtered_Hits.txt", header = TRUE)

#Check for homeostasis proteins in QTL####
#Will be checking using bp and not genetic distance
#We will expect MTF-1 to be listed because it was within our identified candidate
#genes. It will serve as a positive control for these series of filters. 
for (i in 1:nrow(QTL_cM)) {
  #Filter homeo_proteins by QTL's chromosome arm
  x <- homeo_proteins %>% filter(CHROM == QTL_cM$CHROM[i])
  #Filter out proteins that are not within the 3LOD physical range of the QTL
  y <- x %>% filter(LeftEnd_BP > QTL_cM$Left3LOD_phys[i])
  z <- y %>% filter(RightEnd_BP < QTL_cM$Right3LOD_phys[i])
  print(z)
}

# [1] Gene              Class             CHROM             LeftEnd_BP        RightEnd_BP       Chromosome_Recomb Pos_cM           
# <0 rows> (or 0-length row.names)
# [1] Gene              Class             CHROM             LeftEnd_BP        RightEnd_BP       Chromosome_Recomb Pos_cM           
# <0 rows> (or 0-length row.names)
# Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 1   foi   ZIP chr3L    8525244     8534508                 3     26
# 2 MTF-1    TF chr3L    9430312     9439490                 3     29
# [1] Gene              Class             CHROM             LeftEnd_BP        RightEnd_BP       Chromosome_Recomb Pos_cM           
# <0 rows> (or 0-length row.names)
# Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 1 MtnF   Mtn chr3R   24634625    24635165                 3     85
# [1] Gene              Class             CHROM             LeftEnd_BP        RightEnd_BP       Chromosome_Recomb Pos_cM           
# <0 rows> (or 0-length row.names)
# [1] Gene              Class             CHROM             LeftEnd_BP        RightEnd_BP       Chromosome_Recomb Pos_cM           
# <0 rows> (or 0-length row.names)

#As a reminder QTL B on 2L is the first QTL listed in QTL_cM so will be first
#QTL tested
#As expected we do see MTF-1 listed for QTL D.
#We also have foi (ZIP) listed under QTL D and MtnD listed under QTL F. These
#genes were missed when looking for candidate genes

#Compare the location of genes to QTL
#Filter out foi and MtnF
homeo_proteins %>% filter(Gene == "foi" | Gene == "MtnF")
#   Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 1  foi   ZIP chr3L    8525244     8534508                 3     26
# 2 MtnF   Mtn chr3R   24634625    24635165                 3     85

#Filter out QTL D and F
QTL_cM %>% filter(CHROM == "chr3L" | cM == 81.213)
#            CHROM      POS     cM   mlog10p  smooLOD Left3LOD Right3LOD Left3LOD_phys Right3LOD_phys IntSizecM
# chr3L.2541 chr3L  8925893 26.677 13.403906 13.32576   24.327    29.027       8352067        9512583       4.7
# chr3R.3610 chr3R 23906956 81.213  9.948549  9.98031   79.613    84.113      23458263       24716284       4.5

#foi and QTL D: foi is near the left 3 LOD physical end of QTL D. (173,177 bp from end)
#MtnF and QTL F: MtnF is near the right 3 LOD physical end of QTL F. (81,119 bp from end)
#At the time we were identifying candidate genes we had not finalized our QTL
#analysis, which effected where the ends of the peaks were called. This is why
#these genes were not identified when looking for candidates. 

#Check for homeostasis proteins in Filtered out peaks####
for (i in 1:nrow(Filt_Peaks_cM)) {
  x <- homeo_proteins %>% filter(CHROM == Filt_Peaks_cM$CHROM[i])
  y <- x %>% filter(LeftEnd_BP > Filt_Peaks_cM$Left3LOD_phys[i])
  z <- y %>% filter(RightEnd_BP < Filt_Peaks_cM$Right3LOD_phys[i])
  print(z)
}

#foi and MtnF were also identified in the filtered out peaks. This is because they 
#were found in the identified peaks that were filtered out due to being on top
#of another peaks. 

#ZnT77C was identified in two of the filtered out peaks due to the two filtered
#out peaks overlapping. The two peaks it was found in were filtered out because 
#they were in the centromere region

#MtnA and ZnT86D were identified in two peaks that are back to back together. 
#These two peaks were filtered out because they were in the centromere region.

#MtnB-E were identified in the filtered out peak that was next to QTL E which is
#very large. Was thought it could be a remnant of QTL E. 

#Supplemental figure Homeostasis Proteins####
#Create a figure that has the location of each homeostasis protein marked, as well
#as our QTL peaks and filtered out peaks

#Due to how close MtnB-E and Zip42C.1 and Zip42C.2 will be combining into one 
#marker for the plot

Mtn <- homeo_proteins[19,]
Mtn$Gene <- c("MtnB-E")
Mtn$RightEnd_BP <- c(20535301)

Zip42 <- homeo_proteins[10,]
Zip42$Gene <- c("Zip42C.1 & .2")
Zip42$RightEnd_BP <- c(6719533)

#Remove MtnB-E and Zip42C.1 and Zip42C.2 from table of homeostasis proteins
homeo_proteins_combined <- homeo_proteins %>% filter(Pos_cM !=68 & Pos_cM != 69)
homeo_proteins_combined <- homeo_proteins_combined %>% filter(Gene != "Zip42C.1" & Gene != "Zip42C.2")
#Add Mtn and Zip42 to homeo_proteins_combined
homeo_proteins_combined <- rbind(homeo_proteins_combined, Mtn)
homeo_proteins_combined <- rbind(homeo_proteins_combined, Zip42)

#In this graph we will be coloring proteins based upon their class. 
#Not every protein arm has each class so to keep the colors consistent creating
#a "protein" of each class for every arm. 
chr2L_TF <- homeo_proteins_combined[1,]
chr2L_TF$Class <- c("TF")
chr2L_TF$LeftEnd_BP <- NA
chr2L_TF$RightEnd_BP <- NA
chr2L_Mtn <- homeo_proteins_combined[1,]
chr2L_Mtn$Class <- c("Mtn")
chr2L_Mtn$LeftEnd_BP <- NA
chr2L_Mtn$RightEnd_BP <- NA

chr2R_TF <- homeo_proteins_combined[3,]
chr2R_TF$Class <- c("TF")
chr2R_TF$LeftEnd_BP <- NA
chr2R_TF$RightEnd_BP <- NA
chr2R_Mtn <- homeo_proteins_combined[3,]
chr2R_Mtn$Class <- c("Mtn")
chr2R_Mtn$LeftEnd_BP <- NA
chr2R_Mtn$RightEnd_BP <- NA

chr3L_Mtn <- homeo_proteins_combined[6,]
chr3L_Mtn$Class <- c("Mtn")
chr3L_Mtn$LeftEnd_BP <- NA
chr3L_Mtn$RightEnd_BP <- NA

chr3R_TF <- homeo_proteins_combined[7,]
chr3R_TF$Class <- c("TF")
chr3R_TF$LeftEnd_BP <- NA
chr3R_TF$RightEnd_BP <- NA

#Combine the class "proteins" to homeo_proteins_combined
homeo_proteins_combined <- rbind(homeo_proteins_combined, chr2L_TF, chr2L_Mtn,
                                 chr2R_Mtn, chr2R_TF, chr3L_Mtn, chr3R_TF)
head(homeo_proteins_combined)
#     Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 1 ZnT33D   ZnT chr2L   12164241    12166704                 2     46
# 2 ZnT35C   ZnT chr2L   15228728    15247714                 2     51
# 3 ZnT41F   ZnT chr2R    5774776     5782858                 2     55
# 4 ZnT49B   ZnT chr2R   12652805    12655743                 2     67
# 5 ZnT63C   ZnT chr3L    3303708     3314009                 3      7
# 6 ZnT77C   ZnT chr3L   20484293    20494921                 3     47

tail(homeo_proteins_combined)
#       Gene Class CHROM LeftEnd_BP RightEnd_BP Chromosome_Recomb Pos_cM
# 110 ZnT33D    TF chr2L         NA          NA                 2     46
# 111 ZnT33D   Mtn chr2L         NA          NA                 2     46
# 31  ZnT41F   Mtn chr2R         NA          NA                 2     55
# 32  ZnT41F    TF chr2R         NA          NA                 2     55
# 61  ZnT77C   Mtn chr3L         NA          NA                 3     47
# 71  ZnT86D    TF chr3R         NA          NA                 3     50
#Because these will not be plotted their gene name does not matter. 

#Remove some of the filtered out peaks due to overlap
Filthits_to_keep <- c(1,2,4,6,7,8)
#Removed row 3 becuase it overlaps with QTL D.
#Removed row 5 because it is encompassed by row 4's peak.
#Removed row 9 because it overlaps with QTL F. 
nohit <- rep(0,nrow(Filt_Peaks_cM))
nohit[Filthits_to_keep] <- 1
Filter_keep <- Filt_Peaks_cM[nohit==1,]

chr <- c("chr2L", "chr2R", "chr3L", "chr3R")

Homeo_proteins_supplement <- function(GENOMIC_SMOOTH, PROTEINS,PEAKS, NO_PEAKS){
  pdf("Homeostasis Proteins.pdf", width = 20, height = 10)
  for (i in chr) {
    filt_smooth <- GENOMIC_SMOOTH %>% filter(CHROM == i)
    filt_proteins <- PROTEINS %>% filter(CHROM == i)
    filt_QTL <- PEAKS %>% filter(CHROM == i)
    filt_nonhits <- NO_PEAKS %>% filter(CHROM == i)
    plot <- ggplot(filt_smooth)+
      #Mark identified QTL in wheat color rectangle
      geom_rect(data = filt_QTL, mapping = aes(xmin=Left3LOD_phys, xmax=Right3LOD_phys, ymin=-Inf, ymax=40), fill="wheat")+
      #Mark filtered out peaks in gray color rectangle
      geom_rect(data = filt_nonhits, mapping = aes(xmin=Left3LOD_phys, xmax=Right3LOD_phys, ymin=-Inf, ymax=40), fill="gray90", alpha=.5)+
      #Plot smooLOD scores in line
      geom_line(aes(x=POS, y=smooLOD), size=1)+
      #Mark the ends of the genes and color depending on protein class
      geom_rect(data = filt_proteins, mapping = aes(xmin=LeftEnd_BP, xmax=RightEnd_BP, ymin=-Inf, ymax=32, color=Class, fill=Class))+
      #Label each protein rectangle 
      geom_text(data = filt_proteins, aes(x=(RightEnd_BP+LeftEnd_BP)/2, y=35), label=filt_proteins$Gene, angle=90, size=5)+
      #Plot QTL threshold at LOD score 4 in red horizontal line
      geom_hline(yintercept = 4, linetype = 2, colour="red", size=1)+
      xlab("Position (BP)")+
      ylab(expression(-log[10]~italic('(P)')))+
      ggtitle(i)+
      theme_classic()+
      theme(text = element_text(size = 20))
    print(plot)
  }
  graphics.off()
}

Homeo_proteins_supplement(Genomic_smoothed, homeo_proteins_combined, QTL_cM, Filter_keep)
#You will get warning messages saying rows were removed. These are the "proteins" 
#being removed because they do not have a LeftEND_BP or RightEnd_BP

#Sometimes when running this function the page that should have chr2L is blank.
#If you rerun the function it usually will fix this. 

#Catsup will be not be on the black LOD line because our windows did not go that
#far. Catsup is in the centromere region and is likely part of that filtered
#peak. 