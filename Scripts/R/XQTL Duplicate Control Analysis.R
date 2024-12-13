#Introduction####
#Title: Duplicate Control Analysis
#Purpose: Compare duplicate control pools of replicates 5, 6, 8, 11, and 12.
#Check pools have same haplotype frequencies and don't map QTL
#Created: 12/21/22
#Last Edited: 6/22/23
#Packages Needed: 
library(tidyverse) #
library(data.table)#

#Read in Genetic Scan####
# READ in haplotype calls based on GENETIC POSITION
haplo_filename = "allhaps.zinc.015c.txt.gz"
haplo_raw = fread(haplo_filename,header=TRUE)

#You will not need to make any additional directories as long as 
#"allhaps.zinc.015.txt.gz" is in your working directory. 
#This code is very similar to Zinc Analysis R, except for treating replicate pool
#as treatment. Further detail about these codes can be found in that script.
#You do not need to have run Zinc Analysis R to use this script. However, at line
#284, you may load the QTL hit table created in Zinc Analysis R. This is not
#necessary to complete function of this script, and is only meant for easy comparison. 

# TRANSFORM haplotype frequencies
# Same pipe from Zinc Analysis R starting at line 59, but stops before
# grouping and summarizing at line 88. 

haplo_freqs = haplo_raw %>%
  
  # ADD a TRT column that is either C or Z
  mutate(TRT = str_sub(pool, 1, 1)) %>%

  # ADD a REP column that is 01-12
  mutate(REP = str_sub(pool, 2, 3)) %>%

  # ADD a repREP column that is A or B when
  #   there is a Control replicate (e.g. C05A,C05B)
  #   or is an empty cell (when no replicate pools)
  mutate(repREP = str_sub(pool, 4, 5)) %>%
  
  # SWITCH empty cells in repREP to A (effectively
  #   saying for a sample with no replicate pools, its the
  #   A replicate)
  mutate(repREP = recode(na_if(repREP,""), .missing="A")) %>%
  
  # REMOVE the pool column
  select(-c(pool)) 

head(haplo_freqs)
#     chr     pos    cM founder   freq TRT REP repREP
# 1: chrX 1602150 0.924      A1 0.0267   C  01      A
# 2: chrX 1602150 0.924      A2 0.0031   C  01      A
# 3: chrX 1602150 0.924      A3 0.0783   C  01      A
# 4: chrX 1602150 0.924      A4 0.5572   C  01      A
# 5: chrX 1602150 0.924      A5 0.1284   C  01      A
# 6: chrX 1602150 0.924      A6 0.0330   C  01      A

#Isolate Control Duplicates and Format repREP as treatment####
#Filter out zinc treatment and replicates without duplicate control pools
haplo_freqs <- haplo_freqs %>% 
  filter(TRT == "C") %>% 
  filter(REP == "05" | REP == "06" | REP == "08" | REP == "11" |REP == "12")

#Change repREP column into the new TRT column.
#During this analysis we are going treat the two pools as different treatments. 
haplo_freqs <- haplo_freqs %>% select(-c(TRT))
colnames(haplo_freqs)[7] <- c("TRT")
dim(haplo_freqs)
#[1] 422400 7
haplo_freqs
#           chr      pos      cM founder   freq REP TRT
#      1:  chrX  1602150   0.924      A1 0.0281  05   A
#      2:  chrX  1602150   0.924      A2 0.0011  05   A
#      3:  chrX  1602150   0.924      A3 0.0909  05   A
#      4:  chrX  1602150   0.924      A4 0.5335  05   A
#      5:  chrX  1602150   0.924      A5 0.1462  05   A
# ---                                              
# 422396: chr3R 31020051 100.463      A4 0.6303  12   B
# 422397: chr3R 31020051 100.463      A5 0.1329  12   B
# 422398: chr3R 31020051 100.463      A6 0.0385  12   B
# 422399: chr3R 31020051 100.463      A7 0.0647  12   B
# 422400: chr3R 31020051 100.463     AB8 0.0094  12   B

#Summarize/group and arcsin sqrt the frequencies  
haplo_freqs <- haplo_freqs %>% group_by(chr,pos,founder,TRT,REP) %>% 
  summarize(asf=asin(sqrt(mean(freq))),cM=mean(cM))

head(haplo_freqs)
#   chr      pos founder TRT   REP     asf    cM
#   <chr>  <int> <chr>   <chr> <chr> <dbl> <dbl>
# 1 chr2L 590865 A1      A     05    0.139  2.61
# 2 chr2L 590865 A1      A     06    0.136  2.61
# 3 chr2L 590865 A1      A     08    0.149  2.61
# 4 chr2L 590865 A1      A     11    0.150  2.61
# 5 chr2L 590865 A1      A     12    0.166  2.61
# 6 chr2L 590865 A1      B     05    0.143  2.61

#Now haplo_freqs is in the same format it would have been after going through
#line 95 of Zinc Analysis R, except now TRT column is either A or B. 

#Continue with ANOVA Scan in genetic distance analysis####
#Run remaining of the ANOVA Scan with genetic intervals (lines 116-165)
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

# CONVERT the data from chr/pos/founder/TRT/REP/asf columns to
# chr/pos/TRT/REP/A1/A2/A3/A4/A5/A6/A7/AB8 format
# ADD column names 
haplo_freqs_wide = haplo_freqs %>% pivot_wider(names_from=founder,values_from=asf)
colnames(haplo_freqs_wide) = c("CHROM", "POS", "TRT", "REP", "cM", paste0("HAP",1:8))

head(haplo_freqs_wide)
#   CHROM    POS TRT   REP      cM  HAP1  HAP2  HAP3  HAP4  HAP5  HAP6  HAP7  HAP8
#   <chr>  <int> <chr> <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 chr2L 590865 A     05     2.61 0.139 0.143 0.470 0.354 0.252 0.419 0.446 0.489
# 2 chr2L 590865 A     06     2.61 0.136 0.152 0.379 0.379 0.253 0.422 0.485 0.508
# 3 chr2L 590865 A     08     2.61 0.149 0.149 0.387 0.399 0.240 0.432 0.450 0.513
# 4 chr2L 590865 A     11     2.61 0.150 0.161 0.386 0.353 0.277 0.424 0.431 0.549
# 5 chr2L 590865 A     12     2.61 0.166 0.181 0.374 0.328 0.287 0.466 0.415 0.534
# 6 chr2L 590865 B     05     2.61 0.143 0.106 0.380 0.391 0.252 0.407 0.474 0.533

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
# 1 chr2L 590865  2.61 0.00333
# 2 chr2L 617342  2.66 0.00272
# 3 chr2L 643762  2.71 0.00205
# 4 chr2L 670120  2.76 0.00155
# 5 chr2L 696411  2.81 0.00142
# 6 chr2L 722630  2.86 0.00196

write.table(haplo_ANOVA_freq_test_output,"ANOVA_haplotype_freq_scan_cM_Dup_Pool.txt")

#Supplementary Figure 1 Visualize with Manhattan Plot####
make.Manhattan = function(df,Y,myxlab,myylab,mytitle,threshold,ylimit,physical) {
  
  chrlab=c("X","2L","2R","3L","3R")
  myY=sym(Y)
  if(physical==TRUE){
    myX="BPcum"
    totlen = df %>%
      mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
      group_by(Ichr) %>% 
      summarise(chr_len=max(POS)) %>% 
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len)
    
    df2 = df %>% mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
      left_join(totlen, ., by=c("Ichr"="Ichr")) %>%
      arrange(Ichr, POS) %>%
      mutate(BPcum=POS+tot)
    
    temp = df2 %>%
      group_by(Ichr) %>% 
      summarize(center=(max(BPcum) + min(BPcum))/2)
    
    mycenter = as.numeric(temp$center)
    
    ggplot(df2, aes_string(x=myX, y=myY)) +
      ylab(myylab) +
      xlab(myxlab) +
      ggtitle(mytitle) + 
      theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
      # Show all points
      geom_point( aes(color=as.factor(Ichr)), alpha=0.8, size=0.3) +
      scale_color_manual(values = c("grey30", "grey70", "grey30", "grey70", "grey30")) +
      # SNPs of interest
      #		{if(SNPsOfInterest != NULL) geom_point(data=subset(df, ID %in% SNPsOfInterest), color="orange", size=0.8)} +
      # threshold
      {if(threshold != 0) geom_hline(yintercept = threshold, linetype = "dashed", colour = "blue")} +  
      # custom X axis:
      scale_x_continuous(label = chrlab, breaks= mycenter ) +
      scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
      # Custom the theme:
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) +
      theme(panel.border = element_rect(fill = NA, color = "black"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
      theme(legend.position = "none")
    
  }
  else {
    myX="cMcum"
    totlen = df %>%
      mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
      mutate(IIchr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=2,'chr3L'=3,'chr3R'=3)) %>%
      group_by(IIchr) %>% 
      summarise(chr_len=max(cM)) %>% 
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len)
    
    df2 = df %>% mutate(Ichr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=3,'chr3L'=4,'chr3R'=5)) %>%
      mutate(IIchr=recode(CHROM,'chrX'=1,'chr2L'=2,'chr2R'=2,'chr3L'=3,'chr3R'=3)) %>%
      left_join(totlen, ., by=c("IIchr"="IIchr")) %>%
      arrange(Ichr, cM) %>%
      mutate(cMcum=cM+tot)
    
    temp = df2 %>%
      group_by(Ichr) %>% 
      summarize(center=(max(cMcum) + min(cMcum))/2)
    
    mycenter = as.numeric(temp$center)
    
    ggplot(df2, aes_string(x=myX, y=myY)) +
      ylab(myylab) +
      xlab(myxlab) +
      ggtitle(mytitle) + 
      theme(plot.title = element_text(vjust = - 10, hjust=0.025, size=10)) +
      # Show all points
      geom_point( aes(color=as.factor(Ichr)), alpha=0.8, size=0.3) +
      scale_color_manual(values = c("grey30", "grey70", "grey30", "grey70", "grey30")) +
      # SNPs of interest
      #		{if(SNPsOfInterest != NULL) geom_point(data=subset(df, ID %in% SNPsOfInterest), color="orange", size=0.8)} +
      # threshold
      {if(threshold != 0) geom_hline(yintercept = threshold, linetype = "dashed", colour = "blue")} +  
      # custom X axis:
      scale_x_continuous(label = chrlab, breaks= mycenter ) +
      scale_y_continuous(expand = c(0, 0), limits=c(0,ylimit) ) +     # remove space between plot area and x axis
      # Custom the theme:
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank()) + 
      theme(panel.border = element_rect(fill = NA, color = "black"),axis.text=element_text(size=8),axis.title=element_text(size=10)) +
      theme(legend.position = "none")
  }
}
ANOVA_hap_freq_cM = read.table("ANOVA_haplotype_freq_scan_cM_Dup_Pool.txt")
D = make.Manhattan(ANOVA_hap_freq_cM, "mlog10p","location (cM)","-log10(p)","Duplicate Control Pools Haplotype Frequencies", 4,6,TRUE  )
D
#Only 1 window on chromosome X that is above our significant threshold of LOD 4. 
#Check to see if this point is in QTL A
ANOVA_hap_freq_cM %>% arrange(desc(mlog10p))
#      CHROM      POS     cM  mlog10p
# 4841  chrX 13440119 43.174 5.508049 #Not in QTL A which is between 44.7 - 46.0 cM
# 2618 chr3L  9884544 30.527 3.319827 #Not in QTL D which is between 24.3 - 29.0 cM
# 2616 chr3L  9859780 30.427 3.291812
# 2617 chr3L  9872161 30.477 3.269929

#Need to have made this table in Zinc Analysis R.
#This is not essential and only for comparison. 
found_QTL <- read.table("ANOVA_haplotype_freq_scan_cM_QTLhits.txt")
found_QTL
#           CHROM      POS      cM   mlog10p  smooLOD Left3LOD Right3LOD Left3LOD_phys Right3LOD_phys IntSizecM
# chr2L.445  chr2L  7135043  24.809 12.824758 12.50264   23.459    25.659       6829342        7325061      2.20
# chr2R      chr2R 15413541  74.907 11.611377 11.62891   73.807    76.657      15116465       15886141      2.85
# chr3L.2541 chr3L  8925893  26.677 13.403906 13.32576   24.327    29.027       8352067        9512583      4.70
# chr3R.3270 chr3R 18864745  64.213 37.632038 37.52409   63.863    64.513      18717419       18987824      0.65
# chr3R.3610 chr3R 23906956  81.213  9.948549  9.98031   79.613    84.113      23458263       24716284      4.50
# chr3R.3991 chr3R 30848984 100.263 13.994415 13.77635   99.663   100.463      30368845       31020051      0.80
# chrX        chrX 13945053  45.374 17.891269 16.02524   44.724    46.024      13796451       14093253      1.30

#QTL Peak Analysis####
#See Zinc Analysis R lines 181-247 for more detail
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

INFILE = read.table("ANOVA_haplotype_freq_scan_cM_Dup_Pool.txt")
INFILE$CHROM <- as.factor(INFILE$CHROM)
chr_ids <- c("chrX","chr2L","chr2R","chr3L","chr3R")
chr_loess_spans <- matrix(c(0.016, 0.02, 0.02, 0.022, 0.0195),ncol=1)
rownames(chr_loess_spans) <- chr_ids
colnames(chr_loess_spans) <- "loess_span_val"

ll = list()
outINFILE <- NULL
for(chr in levels(INFILE$CHROM)){
  temp = INFILE[INFILE$CHROM==chr,] # EXTRACT 1 chromosome arm of data
  
  fit1 <- loess(temp$mlog10p ~ temp$cM, degree=2, span = chr_loess_spans[rownames(chr_loess_spans)==chr,], family="symmetric")
  temp$smooLOD = fit1$fitted
  
  outINFILE <- rbind(outINFILE,temp)
  
  temphits = temp[find_peaks(temp$smooLOD,20),]
  temphits = temphits[temphits$smooLOD > 7,]
  temphits$Left3LOD = rep(NA,nrow(temphits))
  temphits$Right3LOD = rep(NA,nrow(temphits))
  ll[[chr]] = data.frame(temphits)
}
ll
# $chr2L
# [1] CHROM     POS       cM        mlog10p   smooLOD   Left3LOD  Right3LOD
# <0 rows> (or 0-length row.names)
# 
# $chr2R
# [1] CHROM     POS       cM        mlog10p   smooLOD   Left3LOD  Right3LOD
# <0 rows> (or 0-length row.names)
# 
# $chr3L
# [1] CHROM     POS       cM        mlog10p   smooLOD   Left3LOD  Right3LOD
# <0 rows> (or 0-length row.names)
# 
# $chr3R
# [1] CHROM     POS       cM        mlog10p   smooLOD   Left3LOD  Right3LOD
# <0 rows> (or 0-length row.names)
# 
# $chrX
# [1] CHROM     POS       cM        mlog10p   smooLOD   Left3LOD  Right3LOD
# <0 rows> (or 0-length row.names)

#No peaks found which is supported with the Manhattan plot. 