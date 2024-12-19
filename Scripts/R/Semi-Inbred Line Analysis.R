#Introduction####
#Title: Semi-Inbred Line Analysis on Heavy Metals
#Purpose: The purpose of this script is look at developmental differences
#between the semi-inbred lines derived from Replicate R's zinc selected (9Z) and
#non-selected control (9C) populations. Lines were tested on 10mM ZnCl2 (pilot
#experiments showed that 25mM ZnCl2 was to toxic), H2O, 2mM CuSO4 and 0.2mM CdCl2.
#Will also be creating the figures that make up supplementary material figure 2. 
#Packages 
#Created: 1/13/2023
#Laste Edited: 7/14/2023
#Packages:
library(tidyverse)
library(Rmisc)
library(gridExtra)

#You will need "Semi-Inbred Line Counts.txt"

#As a note: For our analysis we are not focused on individual lines. Instead the lines
#collectively represent either the zinc selected or non selected population.

#Read in Data####
Counts_long <- read.table("Semi-Inbred Line Counts.txt", header = TRUE) 
head(Counts_long)
#   Line Population Treatment      Setup Sex Embryo Total_Flies       Dates Counts
# 1    2         9C      10mM 2021-12-01   M     25           1 X2021.12.12      0
# 2    2         9C      10mM 2021-12-01   M     25           1 X2021.12.13      0
# 3    2         9C      10mM 2021-12-01   M     25           1 X2021.12.14      0
# 4    2         9C      10mM 2021-12-01   M     25           1 X2021.12.15      0
# 5    2         9C      10mM 2021-12-01   M     25           1 X2021.12.16      1
# 6    2         9C      10mM 2021-12-01   M     25           1 X2021.12.17      0

#Line: unique numeric identifier for each inbred line. 9C: 1-22. 9Z: 31-53
#Population: Population that the line was derived from. 9C for Replicate 9's
#control non-selected population and 9Z for Replicate 9's zinc selected population.
#Treatment: Type of media that eggs were put onto. Media was saturated with either
#10mM ZnCl2 (10mM), H2O, 2mM CuSO4 (2mMCuSO4), or 0.2mM CdCl2 (.2mMCdCl2)
#Setup: Date that eggs were collected and put into vials.
#Sex: Sex M or F or flies 
#Embryo: Total number of eggs puts in vials divide by 2 (using the assuming that 
#half of the eggs will M or F). 
#Total_Flies: Total number of flies of that sex that emerged from that vial.
#Counts: number of flies of that sex that emerged from that vial that day. 

#Development Time of Emergence Table####
#Edited "Dates" column into the correct format to get the difference in days. 
#remove"X" at the beginning of the date
Counts_long$Dates <- gsub("X", "",Counts_long$Dates)
#sub out . with -
Counts_long$Dates <- gsub("\\.", "-",Counts_long$Dates)
head(Counts_long)
#   Line Population Treatment      Setup Sex Embryo Total_Flies      Dates Counts
# 1    2         9C      10mM 2021-12-01   M     25           1 2021-12-12      0
# 2    2         9C      10mM 2021-12-01   M     25           1 2021-12-13      0
# 3    2         9C      10mM 2021-12-01   M     25           1 2021-12-14      0
# 4    2         9C      10mM 2021-12-01   M     25           1 2021-12-15      0
# 5    2         9C      10mM 2021-12-01   M     25           1 2021-12-16      1
# 6    2         9C      10mM 2021-12-01   M     25           1 2021-12-17      0

#Get the difference in time in days between the setup and the count date
Diff_Date <- as.numeric(difftime(Counts_long$Dates, Counts_long$Setup, tz="UTC"))
Counts_long$Diff_Date <- Diff_Date
head(Counts_long)
#   Line Population Treatment      Setup Sex Embryo Total_Flies      Dates Counts Diff_Date
# 1    2         9C      10mM 2021-12-01   M     25           1 2021-12-12      0        11
# 2    2         9C      10mM 2021-12-01   M     25           1 2021-12-13      0        12
# 3    2         9C      10mM 2021-12-01   M     25           1 2021-12-14      0        13
# 4    2         9C      10mM 2021-12-01   M     25           1 2021-12-15      0        14
# 5    2         9C      10mM 2021-12-01   M     25           1 2021-12-16      1        15
# 6    2         9C      10mM 2021-12-01   M     25           1 2021-12-17      0        16

#Create a time of emergence table
#Each row will represent a fly that emerged with the time it took to emerge post seutp

#What code is doing is repeating "Diff_Date" (or any other column), by the number
#that is in the counts column. 
Devl_Rep <- rep(Counts_long$Diff_Date, Counts_long$Counts)
Pop_Rep <- rep(Counts_long$Population, Counts_long$Counts)
Setup_Rep <- rep(Counts_long$Setup, Counts_long$Counts)
sex_Rep <- rep(Counts_long$Sex, Counts_long$Counts)
Line_Rep <- rep(Counts_long$Line, Counts_long$Counts)
#Bind the repeated vectors together into one table
Devl <- cbind(Line_Rep, Pop_Rep)
Devl <- as.data.frame(Devl)
Devl$Treatment <- rep(Counts_long$Treatment, Counts_long$Counts)
Devl <- cbind(Devl,Setup_Rep)
Devl <- cbind(Devl, sex_Rep)
Devl <- cbind(Devl, Devl_Rep)
#Give the column's more informative names
colnames(Devl) <- c("Line", "Population","Treatment", "Setup_Date", "Sex", "Emergence_Day")
head(Devl)
#   Line Population Treatment Setup_Date Sex Emergence_Day
# 1    2         9C      10mM 2021-12-01   M            15
# 2    2         9C      10mM 2021-12-01   F            14
# 3    2         9C      10mM 2021-12-01   F            17
# 4    3         9C      10mM 2021-12-01   M            15
# 5    3         9C      10mM 2021-12-01   F            15
# 6    3         9C      10mM 2021-12-01   F            15

#Development Stats Table####
#Create a stats table for the Population on each treatment that has the mean 
#emergence day, standard deviation, standard error, and the 95% Confidence interval.
#These tables can help understand the results from the upcoming analysis. 

Pop_Devl_Stats <- summarySE(Devl, measurevar = "Emergence_Day", 
                            groupvars = c("Population", "Treatment"), na.rm = TRUE)
Pop_Devl_Stats
#   Population Treatment    N Emergence_Day        sd         se         ci
# 1         9C .2mMCdCl2  230     16.621739 2.0624437 0.13599351 0.26795852
# 2         9C      10mM 1005     15.735323 2.3113258 0.07290850 0.14307050
# 3         9C  2mMCuSO4  323     14.801858 2.4893056 0.13850867 0.27249622
# 4         9C       H2O  599      9.838063 0.6635419 0.02711159 0.05324550
# 5         9Z .2mMCdCl2  258     16.771318 2.0927034 0.13028602 0.25656412
# 6         9Z      10mM 1007     14.159881 2.3705573 0.07470260 0.14659077
# 7         9Z  2mMCuSO4  364     14.824176 2.2604196 0.11847818 0.23298978
# 8         9Z       H2O  615      9.731707 0.8277867 0.03337958 0.06555199

#Create a stats table but include sex as a factor 
Sex_Devl_Stats<- summarySE(Devl, measurevar = "Emergence_Day", 
                           groupvars = c("Population", "Treatment", "Sex"), na.rm = TRUE)
Sex_Devl_Stats
#    Population Treatment Sex   N Emergence_Day        sd         se         ci
# 1          9C .2mMCdCl2   F 102     16.833333 1.9654608 0.19460962 0.38605314
# 2          9C .2mMCdCl2   M 128     16.453125 2.1291356 0.18819078 0.37239559
# 3          9C      10mM   F 505     15.716832 2.5097061 0.11168046 0.21941658
# 4          9C      10mM   M 500     15.754000 2.0943392 0.09366169 0.18401988
# 5          9C  2mMCuSO4   F 164     15.365854 2.5015931 0.19534160 0.38572634
# 6          9C  2mMCuSO4   M 159     14.220126 2.3456067 0.18601871 0.36740408
# 7          9C       H2O   F 293      9.744027 0.7759264 0.04533011 0.08921515
# 8          9C       H2O   M 306      9.928105 0.5198307 0.02971674 0.05847578
# 9          9Z .2mMCdCl2   F 109     16.752294 2.0328595 0.19471263 0.38595421
# 10         9Z .2mMCdCl2   M 149     16.785235 2.1421250 0.17548972 0.34678918
# 11         9Z      10mM   F 507     14.071006 2.5333521 0.11251012 0.22104450
# 12         9Z      10mM   M 500     14.250000 2.1920562 0.09803173 0.19260582
# 13         9Z  2mMCuSO4   F 157     15.000000 2.1661735 0.17287947 0.34148665
# 14         9Z  2mMCuSO4   M 207     14.690821 2.3256824 0.16164610 0.31869284
# 15         9Z       H2O   F 323      9.637771 0.8996535 0.05005806 0.09848215
# 16         9Z       H2O   M 292      9.835616 0.7277267 0.04258699 0.08381756

#Development time analysis on all treatments####
#9C and 9Z development time on H2O media
W9C <- Devl %>% filter(Treatment == "H2O" & Population == "9C")
W9Z <- Devl %>% filter(Treatment == "H2O" & Population == "9Z")
t.test(W9C$Emergence_Day, W9Z$Emergence_Day)
# data:  W9C$Emergence_Day and W9Z$Emergence_Day
# t = 2.4732, df = 1169, p-value = 0.01353
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 0.02198497 0.19072728
# sample estimates:
# mean of x mean of y 
# 9.838063  9.731707 

#The 9Z population develops faster then the 9C population. However, the difference
#in means is very small

#9C and 9Z development on 10mM ZnCl2
#Because there is differences in development time on H2O. Will be focused on the
#interaction between treatment and population. 
Zinc_aov_devl <- Devl %>% filter(Treatment == "H2O" | Treatment == "10mM")
Zinc_aov_devl <- aov(lm(Emergence_Day~Population*Treatment*Sex, data = Zinc_aov_devl))
anova(Zinc_aov_devl)
#Response: Emergence_Day
#                           Df  Sum Sq Mean Sq   F value  Pr(>F)    
#Population                  1   892.8   892.8  245.9469 < 2e-16 ***
#Treatment                   1 20130.8 20130.8 5545.6498 < 2e-16 ***
#Sex                         1    13.0    13.0    3.5927 0.05812 .  
#Population:Treatment        1   411.1   411.1  113.2451 < 2e-16 ***
#Population:Sex              1     1.8     1.8    0.4917 0.48322    
#Treatment:Sex               1     1.3     1.3    0.3535 0.55221    
#Population:Treatment:Sex    1     0.8     0.8    0.2136 0.64403    
#Residuals                3218 11681.4     3.6      
#There is a significant difference in the Population*Treatment interaction.
#The 9C population experiences less of a developmental delay due to zinc. 

#9C and 9Z Development on 2mM CuSO4
Copper_aov_devl <- Devl %>% filter(Treatment == "H2O" | Treatment == "2mMCuSO4")
Copper_aov_devl <- aov(lm(Emergence_Day~Population*Treatment*Sex, data = Copper_aov_devl))
anova(Copper_aov_devl)
#Response: Emergence_Day
#Df  Sum Sq Mean Sq   F value    Pr(>F)    
#Population                  1     1.1     1.1    0.4721   0.49209    
#Treatment                   1 11098.6 11098.6 4765.6870 < 2.2e-16 ***
#Sex                         1     7.7     7.7    3.3261   0.06834 .  
#Population:Treatment        1     2.2     2.2    0.9643   0.32624    
#Population:Sex              1    10.2    10.2    4.3607   0.03691 *  
#Treatment:Sex               1    88.8    88.8   38.1506 7.995e-10 ***
#Population:Treatment:Sex    1    18.4    18.4    7.9015   0.00499 ** 
#Residuals                1893  4408.5     2.3
#There is not a significant differences in the Population*Treatment interaction.
#Copper has the same effect on both populations. 

Cadmium_aov_devl <- Devl %>% filter(Treatment == "H2O" | Treatment == ".2mMCdCl2")
Cadmium_aov_devl <- aov(lm(Emergence_Day~Population*Treatment*Sex, data = Cadmium_aov_devl))
anova(Cadmium_aov_devl)
#Response: Emergence_Day
#Df  Sum Sq Mean Sq    F value    Pr(>F)    
#Population                  1     3.6     3.6     2.2054  0.137712    
#Treatment                   1 16649.0 16649.0 10202.7986 < 2.2e-16 ***
#Sex                         1     3.7     3.7     2.2777  0.131431    
#Population:Treatment        1     5.5     5.5     3.3428  0.067674 .  
#Population:Sex              1     1.5     1.5     0.9468  0.330681    
#Treatment:Sex               1    10.9    10.9     6.6818  0.009823 ** 
#Population:Treatment:Sex    1     3.4     3.4     2.0929  0.148166    
#Residuals                1694  2764.3     1.6 
#There is no significant differences in the Population*Treatment interaction.
#Cadmium has about the same effect on both populations.
#As note when look at the stats of the populations, the 9Z population actually
#has a greater developmental delay of 7.04 days, which 9C population is only 
#delayed by 6.78 days. While not significant, this does trend with what is 
#seen within PP-II. 


#Fig 4B & Sup 4B####
#Population Development time on each treatment. 
#Adjust names of treatments
Devl$Treatment <- factor(Devl$Treatment, levels = c("H2O", "10mM",".2mMCdCl2","2mMCuSO4"))

anova_label <- data.frame(Treatment = "H2O", Emergence_Day = 6)

Development_zinc <- Devl%>%filter(Treatment %in% c("H2O","10mM"))%>% ggplot(aes(x=Treatment, y=Emergence_Day))+
  geom_point(position = position_jitterdodge(jitter.height = 0.3, dodge.width = 1),
             aes(group=Population,shape=Population),
             color="grey40",fill="white",
             size=1,show.legend = FALSE)+
  geom_boxplot(aes(fill=Population,color=Population),
               outlier.shape = NA, position = position_dodge(1),linewidth=1,alpha=.5)+
  scale_y_reverse()+
  scale_shape_manual(values = c(16,16))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  scale_color_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = anova_label, aes(x=Treatment, y=Emergence_Day),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Zinc_aov_devl)[4,5], digits = 25)),
            nudge_x = .5, size=5)+
  ggtitle("B. Semi-Inbred Lines")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("10mM ZnCl"[2])))+
  ylab("Development Time (Days)")+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.1,.25))
Development_zinc

zinc_label <- data.frame(Treatment = "10mM", Emergence_Day = 6)
cd_label <- data.frame(Treatment = ".2mMCdCl2", Emergence_Day = 6)
cu_label <- data.frame(Treatment = "2mMCuSO4", Emergence_Day = 6)

Development_all_metals <- Devl%>% ggplot(aes(x=Treatment, y=Emergence_Day))+
  geom_point(position = position_jitterdodge(jitter.height = 0.3, dodge.width = 1),
             aes(group=Population,shape=Population),
             color="grey40",fill="white",
             size=1,show.legend = FALSE)+
  geom_boxplot(aes(fill=Population,color=Population),
               outlier.shape = NA, position = position_dodge(1),linewidth=1,alpha=.5)+
  scale_y_reverse()+
  scale_shape_manual(values = c(16,16))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  scale_color_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = zinc_label, aes(x=Treatment, y=Emergence_Day),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Zinc_aov_devl)[4,5], digits = 25)),
             size=5)+
  geom_text(data = cd_label, aes(x=Treatment, y=Emergence_Day),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Cadmium_aov_devl)[4,5], digits = 2)),
             size=5)+
  geom_text(data = cu_label, aes(x=Treatment, y=Emergence_Day),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Copper_aov_devl)[4,5], digits = 2)),
             size=5)+
  ggtitle("B. Semi-Inbred Lines")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2]),
                            expression("0.2mM CdCl"[2]),expression("2mM CuSO"[4])))+
  ylab("Development Time (Days)")+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.1,.25))
Development_all_metals

#Emergence Tables####
#Get a line for each vial tested and sex with the Total number of flies that emerged. 
Emergence <- Counts_long[-c(8,9,10)]
Emergence <- distinct(Emergence)
head(Emergence)
#    Line Population Treatment      Setup Sex Embryo Total_Flies
# 1     2         9C      10mM 2021-12-01   M   25.0           1
# 18    2         9C      10mM 2021-12-01   F   25.0           2
# 35    3         9C      10mM 2021-12-01   M   17.5           1
# 52    3         9C      10mM 2021-12-01   F   17.5           4
# 69    5         9C      10mM 2021-12-01   M   25.0           6
# 86    5         9C      10mM 2021-12-01   F   25.0           7

#Sum up each vials total number of flies. 
Vial_Emergence <- aggregate(Total_Flies~Line+Population+Treatment+Setup,Emergence, sum)
#sum up each vial total number of eggs
Embryo <- aggregate(Embryo~Line+Population+Treatment+Setup, Emergence, sum)
#Combine Embryo column into vial emergence table.
Vial_Emergence$Embryo <- Embryo$Embryo
#Divide total number of flies by number of eggs, to get the percent emergecne of each vial
Vial_Emergence$Emergence <- Vial_Emergence$Total_Flies/Vial_Emergence$Embryo
#   Line Population Treatment      Setup Total_Flies Embryo Emergence
# 1    2         9C      10mM 2021-12-01           3     50 0.0600000
# 2    3         9C      10mM 2021-12-01           5     35 0.1428571
# 3    5         9C      10mM 2021-12-01          13     50 0.2600000
# 4    6         9C      10mM 2021-12-01          14     50 0.2800000
# 5    7         9C      10mM 2021-12-01          10     50 0.2000000
# 6    8         9C      10mM 2021-12-01          14     50 0.2800000


#Emergence Sex Analysis####
#Can compare the ratio of males to females each other usig Fisher Test
Pop_Total_Sex <- aggregate(Total_Flies~Population+Treatment+Sex, Emergence, sum)
Pop_Total_Sex
#    Population Treatment Sex Total_Flies
# 1          9C .2mMCdCl2   F         102
# 2          9Z .2mMCdCl2   F         109
# 3          9C      10mM   F         505
# 4          9Z      10mM   F         507
# 5          9C  2mMCuSO4   F         164
# 6          9Z  2mMCuSO4   F         157
# 7          9C       H2O   F         293
# 8          9Z       H2O   F         323
# 9          9C .2mMCdCl2   M         128
# 10         9Z .2mMCdCl2   M         149
# 11         9C      10mM   M         500
# 12         9Z      10mM   M         500
# 13         9C  2mMCuSO4   M         159
# 14         9Z  2mMCuSO4   M         207
# 15         9C       H2O   M         306
# 16         9Z       H2O   M         292

#Emergence Analysis####
#Comparison of two population on H2O
H9C <- Vial_Emergence %>% filter(Treatment == "H2O" & Population == "9C")
H9Z <- Vial_Emergence %>% filter(Treatment == "H2O" & Population == "9Z")
t.test(H9C$Emergence, H9Z$Emergence)
# data:  H9C$Emergence and H9Z$Emergence
# t = 0.67042, df = 42.952, p-value = 0.5062
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.06372564  0.12719072
# sample estimates:
# mean of x mean of y 
# 0.5665152 0.5347826 

#9C and 9Z on 10mM ZnCl2
Zinc_aov_emerg <- Vial_Emergence %>% filter(Treatment == "H2O" | Treatment == "10mM")
Zinc_aov_emerg <- aov(lm(Emergence~Population*Treatment, data = Zinc_aov_emerg))
anova(Zinc_aov_emerg)
#Response: Emergence
#                       Df Sum Sq Mean Sq F value    Pr(>F)    
# Population             1 0.0296 0.02963  1.0682    0.3028    
# Treatment              1 1.9582 1.95817 70.5977 1.588e-14 ***
# Population:Treatment   1 0.0008 0.00079  0.0285    0.8661    
# Residuals            172 4.7708 0.02774      
#No significant differences 

#9C and 9Z on 2mM CuSO4
Copper_aov_emerg <- Vial_Emergence %>% filter(Treatment == "H2O" | Treatment == "2mMCuSO4")
Copper_aov_emerg <- aov(lm(Emergence~Population*Treatment, data = Copper_aov_emerg))
anova(Copper_aov_emerg)
# Response: Emergence
#                      Df  Sum Sq Mean Sq F value   Pr(>F)    
# Population            1 0.00044 0.00044  0.0157   0.9007    
# Treatment             1 1.35015 1.35015 48.0470 7.23e-10 ***
# Population:Treatment  1 0.01677 0.01677  0.5968   0.4419    
# Residuals            86 2.41666 0.02810 
#No significant differences

#9C and 9Z on 0.2mM CdCl2
Cadmium_aov_emerg <- Vial_Emergence %>% filter(Treatment == "H2O" | Treatment == ".2mMCdCl2")
Cadmium_aov_emerg <- aov(lm(Emergence~Population*Treatment, data = Cadmium_aov_emerg))
anova(Cadmium_aov_emerg)
# Response: Emergence
#                      Df  Sum Sq Mean Sq F value    Pr(>F)    
# Population            1 0.00317 0.00317  0.1085    0.7426    
# Treatment             1 2.44805 2.44805 83.7887 2.395e-14 ***
# Population:Treatment  1 0.00887 0.00887  0.3035    0.5831    
# Residuals            86 2.51266 0.02922     
#No significant differences 

#Fig 3B &  Sup 3B####
Vial_Emergence$Treatment <- factor(Vial_Emergence$Treatment, levels = c("H2O", "10mM",".2mMCdCl2","2mMCuSO4"))

emerg_anova_label <- data.frame(Treatment = "H2O", Emergence = 1)

Emergence_Zinc <- Vial_Emergence %>% filter(Treatment %in% c("H2O","10mM"))%>% ggplot(aes(x=Treatment, y=Emergence))+
  geom_boxplot(aes(fill=Population),outlier.shape = NA,linewidth=1)+
  geom_point(position = position_jitterdodge(),
             aes(group=Population,shape=Population),
             color="black",fill="white",
             size=3,show.legend = FALSE)+
  scale_shape_manual(values = c(21,21))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = emerg_anova_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ", round(anova(Zinc_aov_emerg)[3,5],4)),
            nudge_x = .5, size=5)+
  ggtitle("B. Semi-Inbred Lines")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("10mM ZnCl"[2])))+
  scale_y_continuous(breaks = c(0.0, 0.25, 0.50, 0.75, 1.0))+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.1,.15))
Emergence_Zinc

zinc_anova_label <- data.frame(Treatment = "10mM", Emergence = 1)
cu_anova_label <- data.frame(Treatment = "2mMCuSO4", Emergence = 1)
cd_anova_label <- data.frame(Treatment = ".2mMCdCl2", Emergence = 1)

Emergence_All_Metals <- Vial_Emergence %>% ggplot(aes(x=Treatment, y=Emergence))+
  geom_boxplot(aes(fill=Population),outlier.shape = NA,linewidth=1)+
  geom_point(position = position_jitterdodge(),
             aes(group=Population,shape=Population),
             color="black",fill="white",
             size=3,show.legend = FALSE)+
  scale_shape_manual(values = c(21,21))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = zinc_anova_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ", round(anova(Zinc_aov_emerg)[3,5],2)),
            size=5)+
  geom_text(data = cd_anova_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ", round(anova(Cadmium_aov_emerg)[3,5],2)),
            size=5)+
  geom_text(data = cu_anova_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ", round(anova(Copper_aov_emerg)[3,5],2)),
            size=5)+
  ggtitle("B. Semi-Inbred Lines")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("10mM ZnCl"[2])))+
  scale_y_continuous(breaks = c(0.0, 0.25, 0.50, 0.75, 1.0))+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.1,.15))
Emergence_All_Metals

