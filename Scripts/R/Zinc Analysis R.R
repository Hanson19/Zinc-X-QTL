#Introduction####
#Title: Zinc Analysis R Code
#Purpose: This script will identify QTL in genetic distance using 
#"allhaps.zinc.015.txt.gz". Numbers from table 1 come from this script, and this
#script will show how Figures 1 and 2 were created. 
#Created: 7/27/22
#Last Edited: 6/23/23
#Packages Needed: 
library(tidyverse)
library(data.table)
library(R.utils)

#You will not need to make any additional directories as long as 
#"allhaps.zinc.015.txt.gz" is in your working directory. 

#Anova Scan with genetic intervals####

#Will need following function
# DEFINE the ANOVA testing function
ANOVA_test_function = function(df) {
  
  # df input has variables HAP / TRT / REP / asf
  df = data.frame(df)
  
  # drop rare haplotypes
  # roughly those with 2% minor allele frequency
  # - asin(sqrt(0.02)) = 0.1418971
  tt = tapply(df$asf,df$HAP,mean) #get the mean of asf for each haplotype
  tt2 = names(tt)[tt > 0.14] #keep haplotypes whose asf greater than 0.14 
  df3 = df %>% filter(HAP %in% tt2) %>% droplevels() #remove haplotypes which are not in tt2
  
  #Run the anova
  out=anova(lm(asf~TRT + HAP + TRT:HAP, data=df3))
  F=out[3,3]/out[4,3]
  df1=out[3,1]
  df2=out[4,1]
  # Returns -log10(P), with P being the TRT:HAP p-value
  -pf(F,df1,df2,lower.tail=FALSE,log.p=TRUE)/log(10)
}

# READ in haplotype calls based on GENETIC POSITION
haplo_filename = "allhaps.zinc.015c.txt.gz"
haplo_raw = fread(haplo_filename,header=TRUE)
head(haplo_raw)
#     chr     pos    cM pool founder   freq
# 1: chrX 1602150 0.924  C01      A1 0.0267
# 2: chrX 1602150 0.924  C01      A2 0.0031
# 3: chrX 1602150 0.924  C01      A3 0.0783
# 4: chrX 1602150 0.924  C01      A4 0.5572
# 5: chrX 1602150 0.924  C01      A5 0.1284
# 6: chrX 1602150 0.924  C01      A6 0.0330

# TRANSFORM haplotype frequencies
# - arcsin sqrt
# - also "pools" each pair of "identical" Control replicates
#   (e.g., "C05A" and "C05B" are pooled)
# - runs very quickly
# - commented lines below code shows first line after each step. 
haplo_freqs = haplo_raw %>%
  #     chr      pos      cM pool founder   freq
  #1:  chrX  1602150   0.924  C01      A1 0.0267
  
  
  # ADD a TRT column that is either C or Z
  mutate(TRT = str_sub(pool, 1, 1)) %>%
  #    chr     pos    cM pool founder   freq TRT
  #1: chrX 1602150 0.924  C01      A1 0.0267   C
  
  
  # ADD a REP column that is 01-12
  mutate(REP = str_sub(pool, 2, 3)) %>%
  #     chr      pos      cM pool founder   freq TRT REP
  #1:  chrX  1602150   0.924  C01      A1 0.0267   C  01
  
  # ADD a repREP column that is A or B when
  #   there is a Control replicate (e.g. C05A,C05B)
  #   or is an empty cell (when no replicate pools)
  mutate(repREP = str_sub(pool, 4, 5)) %>%
  
  # SWITCH empty cells in repREP to A (effectively
  #   saying for a sample with no replicate pools, its the
  #   A replicate)
  mutate(repREP = recode(na_if(repREP,""), .missing="A")) %>%
  #      chr      pos      cM pool founder   freq TRT REP repREP
  #1:  chrX  1602150   0.924  C01      A1 0.0267   C  01      A
  
  # REMOVE the pool column
  select(-c(pool)) %>%
  
  #Re-sort rows of table by chr, then pos, then founder, treatment and replicate.
  group_by(chr,pos,founder,TRT,REP) %>%

  #summarize averages the frequency values across the A and B
  #repREP pools, then arcin sqrt transforms. 
  summarize(asf=asin(sqrt(mean(freq))),cM=mean(cM))
# - Note that after running states:
#   'summarise()' has grouped output by 'chr', 'pos', 'founder', 'TRT'.
#   You can override using the `.groups` argument.
# Don't need to do anything in regards to this message. 

head(haplo_freqs)
# A tibble: 6 x 7
# Groups:   chr, pos, founder, TRT [1]
#   chr      pos founder TRT   REP     asf    cM
#   <chr>  <int> <chr>   <chr> <chr> <dbl> <dbl>
# 1 chr2L 590865 A1      C     01    0.148  2.61
# 2 chr2L 590865 A1      C     02    0.146  2.61
# 3 chr2L 590865 A1      C     03    0.149  2.61
# 4 chr2L 590865 A1      C     04    0.163  2.61
# 5 chr2L 590865 A1      C     05    0.141  2.61
# 6 chr2L 590865 A1      C     06    0.139  2.61

# CONVERT the data from chr/pos/founder/TRT/REP/asf columns to
#   chr/pos/TRT/REP/A1/A2/A3/A4/A5/A6/A7/AB8 format
# ADD column names 
haplo_freqs_wide = haplo_freqs %>% pivot_wider(names_from=founder,values_from=asf)
colnames(haplo_freqs_wide) = c("CHROM", "POS", "TRT", "REP", "cM", paste0("HAP",1:8))
head(haplo_freqs_wide)
# A tibble: 6 x 13
# Groups:   CHROM, POS, TRT [1]
#   CHROM    POS TRT   REP      cM  HAP1  HAP2  HAP3  HAP4  HAP5  HAP6  HAP7  HAP8
#   <chr>  <int> <chr> <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 chr2L 590865 C     01     2.61 0.148 0.142 0.451 0.372 0.256 0.396 0.452 0.504
# 2 chr2L 590865 C     02     2.61 0.146 0.155 0.428 0.369 0.286 0.414 0.420 0.518
# 3 chr2L 590865 C     03     2.61 0.149 0.144 0.433 0.406 0.300 0.415 0.432 0.465
# 4 chr2L 590865 C     04     2.61 0.163 0.150 0.413 0.353 0.273 0.450 0.459 0.479
# 5 chr2L 590865 C     05     2.61 0.141 0.126 0.427 0.373 0.252 0.413 0.460 0.511
# 6 chr2L 590865 C     06     2.61 0.139 0.136 0.392 0.395 0.261 0.418 0.460 0.512

# EXECUTE ANOVA test for each position
haplo_ANOVA_freq_test_output = haplo_freqs_wide %>%
  
  # Reverses back to long format (e.g., haplo_freqs), but with new column names
  pivot_longer(cols = starts_with("HAP"), names_to = "HAP", values_to = "asf") %>%
  
  # REMOVE positions/samples where haplotype freq estimates are NA
  drop_na() %>%
  
  # MAKE some columns factors
  mutate(HAP = as.factor(HAP)) %>%
  mutate(TRT = as.factor(TRT)) %>%
  mutate(REP = as.factor(REP)) %>%
  
  # RE-SORT table
  # WITHIN each chr/pos "make" a sub-table of the 4 TRT/REP/HAP/asf
  #   columns, and the 192 rows (24 samples * 8 founders). This sub-table
  #   is in a 'column' called 'data'
  group_by(CHROM,POS,cM) %>%    # RETAIN cM information
  nest() %>%
  
  # DO ANOVA test at every chr/pos
  mutate(ANOVA_out = map(data, ~ANOVA_test_function(.))) %>%
  
  # REMOVE the 'data column' (so only retain CHROM, POS, and the ANOVA_out result)
  select(-data) %>%
  
  # GET back to a normal column format, so now
  #   have CHROM, POS, ANOVA_out
  unnest(ANOVA_out)

dim(haplo_ANOVA_freq_test_output)
#[1] 5280 4

colnames(haplo_ANOVA_freq_test_output) = c("CHROM","POS","cM","mlog10p")
head(haplo_ANOVA_freq_test_output)
#   CHROM    POS    cM mlog10p
#   <chr>  <int> <dbl>   <dbl>
# 1 chr2L 590865  2.61    5.98
# 2 chr2L 617342  2.66    5.98
# 3 chr2L 643762  2.71    6.11
# 4 chr2L 670120  2.76    6.26
# 5 chr2L 696411  2.81    6.36
# 6 chr2L 722630  2.86    6.32

#We now have a table that for every window in the genome has the Treatment*Haplotype
#interaction p-value in log10. 

#Save output
write.table(haplo_ANOVA_freq_test_output,"ANOVA_haplotype_freq_scan_cM.txt")

#Find QTL Peaks####

find_peaks = function (x, m) {
  shape = diff(sign(diff(x, na.pad = FALSE)))
  pks = sapply(which(shape < 0), FUN = function(i) {
    z = i - m + 1
    z = ifelse(z > 0, z, 1)
    w = i + m + 1
    w = ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks = unlist(pks)
  pks
}

#Read in table that was made in Anova Scan with genetic intervals section
INFILE = read.table("ANOVA_haplotype_freq_scan_cM.txt")

# MAKE 'CHROM' column into a factor
INFILE$CHROM <- as.factor(INFILE$CHROM)

# SET "span" value for loess smoothing
# - using a chromosome arm-specific value
#   so the # of markers for the fit is always
#   approx the same (at ~20-21)
# - assumes number of rows per chromosome arm is
#   chrX=1285, chr2L=1013, chr2R=1009, chr3L=921, chr3R=1052
chr_ids <- c("chrX","chr2L","chr2R","chr3L","chr3R")
chr_loess_spans <- matrix(c(0.016, 0.02, 0.02, 0.022, 0.0195),ncol=1)
rownames(chr_loess_spans) <- chr_ids
colnames(chr_loess_spans) <- "loess_span_val"
chr_loess_spans
#      loess_span_val
#chrX          0.0160
#chr2L         0.0200
#chr2R         0.0200
#chr3L         0.0220
#chr3R         0.0195

# LOOP through chromosomes
# FIND hits
# - Get a per-chromosome arm list of hits (will be in "ll")
# - Also produce a smoothed version of the mlog10p column
#   and append to the input data file as "smooLOD" in
#   "outINFILE"
ll = list()
outINFILE <- NULL
for(chr in levels(INFILE$CHROM)){
  temp = INFILE[INFILE$CHROM==chr,] # EXTRACT 1 chromosome arm of data
  
  # Next line loess smoothing 
  # - degree = 2 (degree of polynomial; can be 1 or 2)
  # - span = a proportion (If # of points is N and span=0.5, then for given
  #   position (x) loess will use the 0.5 * N closest datapoints to x for the fit)
  
  fit1 <- loess(temp$mlog10p ~ temp$cM, degree=2, span = chr_loess_spans[rownames(chr_loess_spans)==chr,], family="symmetric")
  temp$smooLOD = fit1$fitted
  
  # KEEP the newly-smoothed dataset
  outINFILE <- rbind(outINFILE,temp)
  
  temphits = temp[find_peaks(temp$smooLOD,20),]
  temphits = temphits[temphits$smooLOD > 7,]
  temphits$Left3LOD = rep(NA,nrow(temphits))
  temphits$Right3LOD = rep(NA,nrow(temphits))
  ll[[chr]] = data.frame(temphits)
}
head(outINFILE)
#   CHROM    POS    cM  mlog10p  smooLOD
# 1 chr2L 590865 2.609 5.976900 5.940311
# 2 chr2L 617342 2.659 5.978776 6.040660
# 3 chr2L 643762 2.709 6.107760 6.134119
# 4 chr2L 670120 2.759 6.261141 6.220698
# 5 chr2L 696411 2.809 6.358122 6.300563
# 6 chr2L 722630 2.859 6.316684 6.373727
#Our original INFILE, but with additional smoothed LOD score. 

ll[1]
# $chr2L
#      CHROM      POS     cM   mlog10p   smooLOD Left3LOD Right3LOD
# 26   chr2L  1225916  3.859  7.184529  7.025624       NA        NA
# 445  chr2L  7135043 24.809 12.824758 12.502643       NA        NA
# 1000 chr2L 17168849 52.559  8.852928  8.661101       NA        NA
#Pre-list of hits for each arm of chromosome. Can change 1 to any number
#between 1-5 to see the other arms of the chromosomes. 

# WRITE out "outINFILE" 
write.table(outINFILE,"ANOVA_haplotype_freq_scan_cM_smoo.txt")

# LOOP through hits to find 3-LOD CIs
hits = do.call("rbind",ll)
for(i in 1:nrow(hits)){
  chr = hits$CHROM[i]
  temp = outINFILE[outINFILE$CHROM==chr,]
  temp2 = temp %>% filter(cM > hits$cM[i] - 3 & cM < hits$cM[i] + 3 & smooLOD - hits$smooLOD[i] > -3) 
  hits$Left3LOD[i] = min(temp2$cM)
  hits$Right3LOD[i] = max(temp2$cM)
  hits$Left3LOD_phys[i] = min(temp2$POS)
  hits$Right3LOD_phys[i] = max(temp2$POS)
}
hits$IntSizecM = (hits$Right3LOD - hits$Left3LOD)
head(hits)
#            CHROM      POS     cM   mlog10p   smooLOD Left3LOD Right3LOD Left3LOD_phys Right3LOD_phys IntSizecM
# chr2L.26   chr2L  1225916  3.859  7.184529  7.025624    2.609     6.809        590865        2489720      4.20
# chr2L.445  chr2L  7135043 24.809 12.824758 12.502643   23.459    25.659       6829342        7325061      2.20
# chr2L.1000 chr2L 17168849 52.559  8.852928  8.661101   51.709    53.209      16125400       18371012      1.50
# chr2R      chr2R 15413541 74.907 11.611377 11.628915   73.807    76.657      15116465       15886141      2.85
# chr3L.2504 chr3L  8472330 24.827 12.086338 12.331639   24.177    27.777       8316145        9200426      3.60
# chr3L.2541 chr3L  8925893 26.677 13.403906 13.325758   24.327    29.027       8352067        9512583      4.70


#To see the smoothing process run the lines below
par(mar=c(1,1,1,1))#Remove this line if you just want to see each individually 
par(mfrow=c(5,1)) 
for(chr in levels(INFILE$CHROM)){
  temp = outINFILE[outINFILE$CHROM==chr,]
  plot(temp$POS,temp$mlog10p,cex=0.5,main=chr)
  lines(temp$POS,temp$smooLOD,col=2,lwd=1.5)
}

#Our code automatically found 16 hits, we then manually went through and filtered the
#hits. Our manual filtration defines which peaks we will be focusing on for the rest
#of our analysis. This doesn't mean the peaks we ignore are not valid.
#Reasons removal:
#Unsure if real peak: ex. row 2
#Near centromere: ex. row 3,7,8,9,10
#Peak on top of peak: ex. row 5, 14
#Proximity another peak: row 12 (next to huge peak of QTL E, could be remnant of E) 

rows_to_keep <- c(2,4,6,11,13,15,16)
qtl <- rep(0,nrow(hits))
qtl[rows_to_keep] <- 1
only_hits <- hits[qtl==1,]
only_hits
#            CHROM      POS      cM   mlog10p  smooLOD Left3LOD Right3LOD Left3LOD_phys Right3LOD_phys IntSizecM
# chr2L.445  chr2L  7135043  24.809 12.824758 12.50264   23.459    25.659       6829342        7325061      2.20
# chr2R      chr2R 15413541  74.907 11.611377 11.62891   73.807    76.657      15116465       15886141      2.85
# chr3L.2541 chr3L  8925893  26.677 13.403906 13.32576   24.327    29.027       8352067        9512583      4.70
# chr3R.3270 chr3R 18864745  64.213 37.632038 37.52409   63.863    64.513      18717419       18987824      0.65
# chr3R.3610 chr3R 23906956  81.213  9.948549  9.98031   79.613    84.113      23458263       24716284      4.50
# chr3R.3991 chr3R 30848984 100.263 13.994415 13.77635   99.663   100.463      30368845       31020051      0.80
# chrX        chrX 13945053  45.374 17.891269 16.02524   44.724    46.024      13796451       14093253      1.30

#WRITE out hit table
write.table(only_hits,"ANOVA_haplotype_freq_scan_cM_QTLhits.txt")
#This table is where Table 1's information (excluding Genes, and Candidate Genes
#columns) is from

#Optional: Saving Filtered out hits
#This is needed for running XQTL Location Homeostasis Proteins
Filt_hits <- hits[qtl==0,]
write.table(Filt_hits, "ANOVA_haplotype_freq_scan_cM_Filtered_Hits.txt")

#Figure 1 QTL Plot####
ANOVA_hap_freq_cM_smoo = read.table("ANOVA_haplotype_freq_scan_cM_smoo.txt")
#Factor chromosomes so chrX is first one. 
ANOVA_hap_freq_cM_smoo$CHROM <- factor(ANOVA_hap_freq_cM_smoo$CHROM, levels = c("chrX","chr2L","chr2R","chr3L","chr3R"))

#Data Points for Labels
QA <- data.frame(cM=45.374, smooLOD=17.5, CHROM=factor("chrX", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QB <- data.frame(cM=24.809, smooLOD=14.25, CHROM=factor("chr2L", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QC <- data.frame(cM=74.907,smooLOD=13.25, CHROM=factor("chr2R", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QD <- data.frame(cM=26.677, smooLOD=15, CHROM=factor("chr3L", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QE <- data.frame(cM=64.213, smooLOD=39.25, CHROM=factor("chr3R", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QF <- data.frame(cM=81.213, smooLOD=11.5, CHROM=factor("chr3R", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))
QG <- data.frame(cM=99.5, smooLOD=15.5, CHROM=factor("chr3R", levels = c("chrX","chr2L","chr2R","chr3L","chr3R")))

Color <- c("red4", "red", "red4", "red", "red4")

QTL_Plot <- ggplot(ANOVA_hap_freq_cM_smoo)+
  geom_line(aes(x=cM, y=smooLOD, group=CHROM, color=CHROM), linewidth =2)+
  scale_color_manual(values = Color)+
  facet_wrap(~CHROM, scales = "free_x", ncol=5, strip.position = "bottom")+
  geom_text(data = QA, aes(cM, smooLOD), label="A", size=10)+
  geom_text(data = QB, aes(cM, smooLOD), label="B", size=10)+
  geom_text(data = QC, aes(cM, smooLOD), label="C", size=10)+
  geom_text(data = QD, aes(cM, smooLOD), label="D", size=10)+
  geom_text(data = QE, aes(cM, smooLOD), label="E", size=10)+
  geom_text(data = QF, aes(cM, smooLOD), label="F", size=10)+
  geom_text(data = QG, aes(cM, smooLOD), label="G", size=10)+
  geom_hline(yintercept = 4, linetype = 2, colour="black", size=2)+
  xlab("Position (cM)")+
  ylab(expression(-log[10]~italic('(P)')))+
  theme_classic()+ 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.position = "none")
QTL_Plot

ggsave("Figure 1 QTL Plot.png", QTL_Plot, width = 15, height = 6, units = "in")

#
#Figure 2 QTL Haplotype Frequency plots for genetic distance####
#Plot founder haploty frequencies
library(tidyverse)
library(ggpubr)
library(gridExtra)

# SET filenames
# READ in haplotype frequency data in Genetic space
# READ in QTL hit information
input_haplo_freq_name <- "allhaps.zinc.015c.txt.gz"
input_qtl_name <- "ANOVA_haplotype_freq_scan_cM_QTLhits.txt"
LODs <- read.table("ANOVA_haplotype_freq_scan_cM_smoo.txt", header = TRUE)
haplo_freq = read.table(input_haplo_freq_name,header=TRUE)
QTL = read.table(input_qtl_name,header=TRUE)

# GET average haplotype frequencies (per founder per position)
#   for each treatment (Selection / Control)
# THEN convert to a "Selection minus Control" difference
haplo_freq_diff_to_plot = haplo_freq %>%
  
  # GENERATE a TRT (treatment) vector that is simply a single
  #   letter pulled from the replicate name
  # REMOVE 'pool' column
  mutate(TRT = substr(pool,1,1)) %>%
  select(-pool) %>%
  #     chr     pos    cM founder   freq TRT
  #1   chrX 1602150 0.924      A1 0.0267   C
  #2   chrX 1602150 0.924      A2 0.0031   C  
  
  # GROUP
  # CALCULATE MEAN of haplotype frequencies (over replicates)
  #   for each founder and marker interval
  group_by(chr,cM,founder,TRT) %>%
  summarize(freq = mean(freq)) %>%
  # chr      cM founder TRT     freq
  #   <chr> <dbl> <chr>   <chr>  <dbl>
  #1 chr2L  2.61 A1      C     0.0216
  #2 chr2L  2.61 A1      Z     0.0208
  
  # REFORMAT to put the treatment-specific frequencies (per
  #   founder and position) in separate columns
  # CALCULATE DIFFERENCE (Selection minus Control)
  # REMOVE 'C' and 'Z' columns
  pivot_wider(names_from=TRT,values_from=freq) %>%
  mutate(Diff = Z-C) %>%
  select(-c(Z,C))

head(haplo_freq_diff_to_plot)
#   chr      cM founder      Diff
#   <chr> <dbl> <chr>       <dbl>
# 1 chr2L  2.61 A1      -0.000796
# 2 chr2L  2.61 A2       0.00105 
# 3 chr2L  2.61 A3      -0.0192  
# 4 chr2L  2.61 A4       0.000594
# 5 chr2L  2.61 A5       0.0127  
# 6 chr2L  2.61 A6       0.0294 
#End of pipeline we have ended up with differences in frequencies for each founder at each position 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

LODPLOT = list()
FREQPLOT = list()
scaleFUN <- function(x) sprintf("%3.1f", x)
mylabs = c("A","B","C","D","E","F","G")

#Arrange QTL order so X chromosome comes first
QTL$POS <- factor(QTL$POS, levels = c(13945053, 7135043, 15413541,8925893,18864745, 23906956, 30848984))
QTL <- QTL %>% arrange(POS)
QTL
#            CHROM      POS      cM   mlog10p  smooLOD Left3LOD Right3LOD Left3LOD_phys Right3LOD_phys IntSizecM
# chrX        chrX 13945053  45.374 17.891269 16.02524   44.724    46.024      13796451       14093253      1.30
# chr2L.445  chr2L  7135043  24.809 12.824758 12.50264   23.459    25.659       6829342        7325061      2.20
# chr2R      chr2R 15413541  74.907 11.611377 11.62891   73.807    76.657      15116465       15886141      2.85
# chr3L.2541 chr3L  8925893  26.677 13.403906 13.32576   24.327    29.027       8352067        9512583      4.70
# chr3R.3270 chr3R 18864745  64.213 37.632038 37.52409   63.863    64.513      18717419       18987824      0.65
# chr3R.3610 chr3R 23906956  81.213  9.948549  9.98031   79.613    84.113      23458263       24716284      4.50
# chr3R.3991 chr3R 30848984 100.263 13.994415 13.77635   99.663   100.463      30368845       31020051      0.80


Founder_Color <- c("#F0E442", "#555555", "#E69F00", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#CC79A7") 

for(i in 1:7){
  locchr=QTL$CHROM[i]
  locpos = QTL$cM[i]
  #Plotting 5cM up and down from QTL peak
  left = locpos - 5
  right = locpos + 5
  myLODstemp = LODs %>% filter(CHROM==locchr & cM > left & cM < right)
  myFstemp =  haplo_freq_diff_to_plot %>% filter(chr==locchr & cM > left & cM < right)
  #Plots the QTL peak which is the top plot. 
  tempLODPLOT = ggplot(myLODstemp, aes(x=cM, y=smooLOD)) + 
    geom_line(size=1) + #change line size
    ylab(expression(-log[10]~italic('P'))) +
    xlab(paste0("Location (cM) on ", QTL$CHROM[i])) +
    geom_vline(xintercept = locpos, linetype = "dashed", colour = "blue", size=.75) + 
    scale_y_continuous(labels=scaleFUN) + 
    # Custom the theme:
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    theme(axis.text=element_text(size=5),axis.title=element_text(size=10)) +  
    theme(legend.position = "none") + 
    theme(axis.title.x=element_blank()) + 
    theme(axis.text.x=element_blank()) +
    annotate("text", label=mylabs[i], x=-Inf, y=Inf, hjust=0, vjust=1, size=5) #changes QTL label
  
  #Plots the founder haplotype frequency change
  tempFREQPLOT = ggplot(myFstemp, aes(x=cM, y=100*Diff, group=founder)) +
    geom_smooth(data=myFstemp, method = "loess",
                aes(color=founder), se=FALSE, show.legend=TRUE)+
    ylab("Freq Change (%)") +
    xlab(paste0("Location (cM) on ", QTL$CHROM[i])) +
    geom_vline(xintercept = locpos, linetype = "dashed", colour = "blue", size=.75) +  
    scale_y_continuous(labels=scaleFUN) +
    scale_color_manual(values = Founder_Color)+
    # Custom the theme:
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
    theme(axis.text=element_text(size=5),axis.title=element_text(size=10)) +
    theme(legend.position = "none")
  
  LODPLOT[[i]] = tempLODPLOT 
  FREQPLOT[[i]] = tempFREQPLOT
  #Create a legend
  if(i == 1){
    Wantlegend = ggplot(myFstemp, aes(x=cM, y=100*Diff, group=founder)) +
      geom_line(aes(color=founder), size=5) +
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
      theme(axis.text=element_text(size=15),axis.title=element_text(size=17)) +
      theme(legend.direction = "horizontal", legend.text = element_text(size = 10),legend.title = element_text(size=15),
            legend.position = "bottom", legend.box = "verical")+
      scale_color_manual(values = Founder_Color)+
      labs(color = "Founder")
    Mylegend = g_legend(Wantlegend)
  }
  
}

#Combine QTL Lod and frequency plots into one figure
lay = rbind(c(1,2),c(8,9), c(3,4),c(10,11), c(5,6),c(12,13), c(7,15),c(14,15))
FreqPlot <- grid.arrange(LODPLOT[[1]], LODPLOT[[2]], LODPLOT[[3]], LODPLOT[[4]], LODPLOT[[5]], LODPLOT[[6]], LODPLOT[[7]],
                         FREQPLOT[[1]], FREQPLOT[[2]], FREQPLOT[[3]], FREQPLOT[[4]], FREQPLOT[[5]], FREQPLOT[[6]], FREQPLOT[[7]],
                         Mylegend, heights = c(1,2,1,2,1,2,1,2), layout_matrix = lay)
#Save figure
ggsave("Figure 2 Haplotype Freqency Plot Smooth.png", FreqPlot, width = 7, height = 8, units = "in")

