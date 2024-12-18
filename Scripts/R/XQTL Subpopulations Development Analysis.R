#Introduction####
#Title: X-QTL Subpopulation Development Resistance Analysis
#Purpose: Analyze development resistance data collected from generation 2 X-QTL
#Subpopulation flies. Will be analyzing this data both for emergence and development
#time.
#Created: 1/9/2023
#Last Edited: 12/18/2024
#Packages:
library(Rmisc)
library(tidyverse)

#You will need "PP-I G2 Embryo Counts.txt". 

#Note:
#At time that this data was collected, we were in middle of Covid-19 pandemic, 
#and we were not aware of development time being a possible phenotype. 
#Due to this flies were not counted every day.However, when they were counted the 
#control nonselected population was always counted at the same time as zinc 
#selected population under the same treatment. The populations that developed on 
#control media where counted 3 days after the first flies emerged on control 
#media, so the majority of flies for control treatment were collected on the same day.

#Read in data####
G2_all <- read.table("PP-I G2 Embryo Counts.txt", header = TRUE)
names(G2_all)[2] <- "Population" 
head(G2_all)
#   Replicate Population  Treatment Bottle Setup_Date Experimentor Embryos Sex X2021.01.15 X2021.01.16 X2021.01.18
# 1         9       Zinc 25mM_ZnCl2      1 2021-01-03          SJM      25   M           0           0           0
# 2         9       Zinc 25mM_ZnCl2      2 2021-01-03          SJM      25   M           0           0           1
# 3         9       Zinc 25mM_ZnCl2      3 2021-01-03          SJM      25   M           0           0           1
# 4         9       Zinc 25mM_ZnCl2      4 2021-01-03          SJM      25   M           0           0           1
# 5         9       Zinc 25mM_ZnCl2      5 2021-01-03          SJM      25   M           0           0           0
# 6         9       Zinc 25mM_ZnCl2      6 2021-01-03          SJM      25   M           0           0           1
#   X2021.01.19 X2021.01.20 X2021.01.22 X2021.01.24 X2021.01.25 X2021.01.26 X2021.01.27 Total_Flies
# 1          NA           0           1           0           0           0           0           1
# 2          NA           0           1           1           0           0           0           3
# 3          NA           2           1           0           0           0           0           4
# 4          NA           2           0           0           2           0           0           5
# 5          NA           1           0           0           0           0           0           1
# 6          NA           1           0           0           0           0           0           2

#Replicate: What replicate populations, 9 or 10, did flies come from. 
#Population: Zinc - Zinc selected population, Control - Control Non-selected Population
#Treatment: Eggs were put on either control H2O media or 25mM ZnCl2 media. 
#Bottle: unique vial ID for each Replicate/Population/Treatment group. 
#1-10: Replicate 9's Zinc selected Population on 25mM ZnCl2
#1-10: Replicate 9's Control non-selected  Population on 25mM ZnCl2
#1-10: Replicate 10's Zinc selected Population on 25mM ZnCl2
#1-10: Replicate 10's Control non-selected  Population on 25mM ZnCl2
#1-3: Replicate 9's zinc selected population on control H2O media
#1-3: Replicate 9's control non-selected population on control H2O media
#1-3: Replicate 10's zinc selected population on control H2O media
#1-3: Replicate 10's control non-selected population on control H2O media
#Setup_Date: Date eggs were collected and put on media 
#Experimentor: Person who collected eggs
#Embryos: Number of eggs put in each vial. *Because data is split by sex and using
#assumption that M:F eggs 50:50 when collected this column is the number of eggs
#put in the vials divided by two. 
#Sex: Emerged fly was male or female 
# X2021.01.15:X2021.01.27: Date flies were counted. 
#Total_Flies: Total number of flies that emerged 

#Clean Up Data. 
G2_all %>% filter(Total_Flies > 25)
#Bottle 1 of Replicate 9's control non-selected population treated on control H2O
#media had more then 50 flies emerge from the bottle indicating that more then 50
#eggs were put on the media. We have to remove this vial from the analysis because 
#we cannot accurately calculate emergence nor are the eggs in similar environment as 
#the other bottle counterparts given differences in eggs numbers
#The three other vials that had a sex count greater then, did not have more then 
#50 flies emerge when looking at the number of flies of the opposite sex.

G2 <- G2_all %>% filter(Total_Flies < 30)


#Emergence Analysis####
#To do this analysis we look at individual vial emergence. 
#First going to do an anova considering all our variables to see which are ones we do not need to consider. 
Emerg <- aov(lm(Total_Flies~Replicate*Population*Treatment*Experimentor*Sex, data = G2))
anova(Emerg)
#Response: Total_Flies
#                           Df Sum Sq Mean Sq   F value    Pr(>F)    
#Replicate                   1    0.4     0.4    0.0908 0.7639041    
#Population                  1   78.0    78.0   16.0583 0.0001323 ***
#Treatment                   1 6174.5  6174.5 1271.2181 < 2.2e-16 ***
#Experimentor                1    0.0     0.0    0.0011 0.9736095    
#Sex                         1    1.2     1.2    0.2442 0.6224553    
#Replicate:Population        1    6.4     6.4    1.3196 0.2539284    
#Replicate:Treatment         1    7.7     7.7    1.5948 0.2101316    
#Population:Treatment        1   10.1    10.1    2.0696 0.1539781    
#Treatment:Experimentor      1    6.6     6.6    1.3563 0.2474754    
#Replicate:Sex               1    6.2     6.2    1.2839 0.2604009    
#Population:Sex              1    0.4     0.4    0.0797 0.7784433    
#Treatment:Sex               1   75.7    75.7   15.5797 0.0001637 ***
#Experimentor:Sex            1    8.4     8.4    1.7338 0.1915065    
#Replicate:Population:Sex    1    5.6     5.6    1.1629 0.2839447    
#Replicate:Treatment:Sex     1   20.4    20.4    4.1924 0.0437311 *  
#Population:Treatment:Sex    1   41.1    41.1    8.4667 0.0046283 ** 
#Treatment:Experimentor:Sex  1   10.0    10.0    2.0645 0.1544745    
#Residuals                  84  408.0     4.9

#Due to Sex and experimenter not significant factors. Will not consider
#these factors in later analysis

#Merge Male and Females counts together.
G2_merged <- aggregate(Total_Flies~Replicate+Population+Treatment+Bottle+Setup_Date, G2, sum)
unique(G2$Embryos)
#[1] 25
#Divide Total Flies by 50 to create emergence column
G2_merged$Emergence <- G2_merged$Total_Flies/50
head(G2_merged)
#   Replicate Population  Treatment Bottle Setup_Date Total_Flies Emergence
# 1         9    Control 25mM_ZnCl2      1 2021-01-03           4      0.08
# 2         9       Zinc 25mM_ZnCl2      1 2021-01-03           3      0.06
# 3         9       Zinc        H20      1 2021-01-03          49      0.98
# 4         9    Control 25mM_ZnCl2      2 2021-01-03           1      0.02
# 5         9       Zinc 25mM_ZnCl2      2 2021-01-03           8      0.16
# 6         9    Control        H20      2 2021-01-03          40      0.80

#ANOVA comparing within a treatment, due to large differences in number of vials 
#between treatments. 

#Zinc Treatment
G2_zinc <- G2_merged %>% filter(Treatment == "25mM_ZnCl2")
Zinc_Emerg <- aov(lm(Emergence~Replicate*Population, data = G2_zinc))
anova(Zinc_Emerg)
#Response: Emergence
#                     Df Sum Sq Mean Sq F value   Pr(>F)   
#Replicate             1   5.63   5.625  0.5677 0.456072   
#Population            1  87.02  87.025  8.7830 0.005363 **
#Replicate:Population  1   0.62   0.625  0.0631 0.803123   
#Residuals            36 356.70   9.908         
#There is a significant difference between the control non-selected population
#and the zinc selected population on zinc media

#Control H2O Treatment
G2_H2O <- G2_merged %>% filter(Treatment != "25mM_ZnCl2")
Control_Emerg <- aov(lm(Emergence~Replicate*Population, data = G2_H2O))
anova(Control_Emerg)
#Response: Total_Flies
#                     Df   Sum Sq   Mean Sq F value Pr(>F)
#Replicate             1 0.016878 0.0168776  2.4751 0.1597
#Population            1 0.000095 0.0000948  0.0139 0.9094
#Replicate:Population  1 0.007585 0.0075852  1.1124 0.3266
#Residuals             7 0.047733 0.0068190 
#There is not a difference between the two populations on control H2O media. 

#Figure 3A Population Emergence####
#Get the mean and CI for emergence grouped by population and treatment, because
#replicate wasn't significant in above anovas not including it in summary.
G2_merged$Replicate <- paste("R", G2_merged$Replicate, sep = "")
G2_merged$Replicate <- as.factor(G2_merged$Replicate)
G2_merged$Treatment[G2_merged$Treatment == "25mM_ZnCl2"] <- c("25mM ZnCl2")
G2_merged$Treatment <- gsub("0", "O", G2_merged$Treatment)
G2_merged$Treatment <- factor(G2_merged$Treatment, levels = c("H2O","25mM ZnCl2"))
G2_merged$Population[G2_merged$Population == "Control"] <- c("Control Non-Selected")
G2_merged$Population[G2_merged$Population == "Zinc"] <- c("Zinc Selected")

subpop_emerg_Anova <- aov(lm(Emergence~Replicate*Population*Treatment, data = G2_merged))
subpop_emerg_table <- anova(subpop_emerg_Anova)
subpop_emerg_lable <- data.frame(Treatment = "H2O", Emergence = 1.15)

water_trt <- t.test(G2_merged %>% filter(Treatment == "H2O" & Population == "Zinc Selected") %>% pull(Emergence),
                    G2_merged %>% filter(Treatment == "H2O" & Population != "Zinc Selected") %>% pull(Emergence))

zinc_trt <- t.test(G2_merged %>% filter(Treatment == "25mM ZnCl2" & Population == "Zinc Selected") %>% pull(Emergence),
                   G2_merged %>% filter(Treatment == "25mM ZnCl2" & Population != "Zinc Selected") %>% pull(Emergence))

water_label <- data.frame(Treatment = "H2O", Emergence = 1.05)
zinc_label <- data.frame(Treatment = "25mM ZnCl2", Emergence = 1.05)

Subpop_Emerg_Plot <- G2_merged %>% ggplot(aes(x=Treatment, y=Emergence))+
  geom_boxplot(aes(fill=Population), outlier.shape = NA, linewidth=1)+
  geom_point(position = position_jitterdodge(),
             aes(group=Population,shape=Replicate),
             color="black",fill="white",
             size=3)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = subpop_emerg_lable, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ", round(subpop_emerg_table[6,5],4)), 
            nudge_x = .5, size=5)+
  geom_text(data = water_label, aes(x=Treatment, y=Emergence),
            label=paste("T-test, ", round(water_trt$p.value,4)), size=5)+
  geom_text(data = zinc_label, aes(x=Treatment, y=Emergence),
            label=paste("T-test, ", round(zinc_trt$p.value,4)), size=5)+
  ggtitle("A. X-QTL Subpopulations")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2])))+
  scale_y_continuous(breaks = c(0.0, 0.25, 0.50, 0.75, 1))+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.15, .3))

Subpop_Emerg_Plot

#Sex Emergence Analysis####
#For Sex we have to do it by counts and not percentage because we cannot assume 
#that we put 25 male eggs on media and 25 female eggs 

#Below Plots helps with visualization/conclusion from Fisher.test
Stats_TF_Sex <- summarySE(G2, measurevar = "Total_Flies", 
                          groupvars = c("Population", "Treatment","Sex"), 
                          na.rm = TRUE)
Stats_TF_Sex %>% ggplot(aes(x=Population, y=Total_Flies, color=Sex))+
  geom_errorbar(aes(ymin=Total_Flies-ci, ymax=Total_Flies+ci), width=0.3, position = position_dodge(.25))+
  geom_point(size=2, position = position_dodge(.25))+
  facet_wrap(~Treatment, scales = "free_x",  ncol = 5, strip.position = "bottom")+
  ylab("Total Flies")+
  theme_classic()+theme(legend.position = c(.1, .8))

#Fisher.test
#We will compare male and female counts within a population across the two 
#treatments, and within treatments across the two populations.

#Add up the emerged flies for each sex based off of population and treatment.
G2_sex <- aggregate(Total_Flies~Population+Treatment+Sex, G2, sum)
#   Population  Treatment Sex Total_Flies
# 1    Control 25mM_ZnCl2   F          44
# 2       Zinc 25mM_ZnCl2   F          89
# 3    Control        H20   F         106
# 4       Zinc        H20   F         116
# 5    Control 25mM_ZnCl2   M          37
# 6       Zinc 25mM_ZnCl2   M          51
# 7    Control        H20   M         111
# 8       Zinc        H20   M         145

#comparing M to F counts within zinc population between two treatments 
#zinc Population
fisher.test(matrix(c(G2_sex[4,4],G2_sex[2,4],G2_sex[8,4],G2_sex[6,4]), ncol = 2))
#p-value = 0.0003505

#Control Population
fisher.test(matrix(c(G2_sex[3,4], G2_sex[1,4],G2_sex[7,4],G2_sex[5,4]), ncol = 2))
#p-value = 0.4359

#So when comparing within a population the M to F ratio does not change in the 
#Control non-selected population, but it does for the zinc selected population. 
#Looking at the G2 Male and Female counts  for the zinc population, it actually 
#switches which sex more in abundance. On control media there are more males then
#females, but on zinc media there are more females then males. 

#Comparing M to F counts within same treatment but different population. 
#25mM ZnCl2 Treatment
fisher.test(matrix(c(G2_sex[1,4], G2_sex[5,4], G2_sex[2,4], G2_sex[6,4]), ncol = 2))
#p-value = 0.2003

#H2O Treatment
fisher.test(matrix(c(G2_sex[3,4], G2_sex[7,4], G2_sex[4,4],G2_sex[8,4]),ncol = 2))
#p-value = 0.3577
#When comparing between populations within the same treatment there is not 
#a significant difference for either treatment.

#Time of Emergence Table####
#Development time for Subpopulation flies is not as accurate as it is for the other
#phenotyping population analysis due to Covid-19 Pandemic and not being able to 
#to come in every day and count flies. 
#We also didn't realize that development time would be a phenotype until 
#retroactively re-analyzing the data, so we waited until most of the control flies
#had emerged before counting. 
#While this is a short coming of this experiment, we did make sure that every 
#time we counted a treatment condition both populations were counted. 
#Due to not counting every day until more towards the end of the experiment, 
#Replicates 9 and 10 populations can have off counts for example; 
#For a day of counting R9 could be at day 12, but R10 could be at day 11.
#Because of this,we cannot compare them based on development time equally.

#In this section we will be creating a time of emergence table. Each row in this
#table will represent a fly that emerged during the experiment, and will include
#identifying information such as which replicate and population it came from, sex,
#day it emerged and what treatment (H2O or 25mM_ZnCl2) it developed on. 
#This is very similar to the time of death table made in Subpopulation Adult Analysis, 
#except instead of having a line recording each fly's death, it is each fly's 
#emergence day. 
#Using individual flies and their development, rather than using a vial as a whole,
#gives us more power because for our analysis. Instead of working with 102 vials
#we will be working with 699 flies. 

head(G2)
#   Replicate Population  Treatment Bottle Setup_Date Experimentor Embryos Sex X2021.01.15 X2021.01.16
# 1         9       Zinc 25mM_ZnCl2      1 2021-01-03          SJM      25   M           0           0
# 2         9       Zinc 25mM_ZnCl2      2 2021-01-03          SJM      25   M           0           0
# 3         9       Zinc 25mM_ZnCl2      3 2021-01-03          SJM      25   M           0           0
# 4         9       Zinc 25mM_ZnCl2      4 2021-01-03          SJM      25   M           0           0
# 5         9       Zinc 25mM_ZnCl2      5 2021-01-03          SJM      25   M           0           0
# 6         9       Zinc 25mM_ZnCl2      6 2021-01-03          SJM      25   M           0           0
#   X2021.01.18 X2021.01.19 X2021.01.20 X2021.01.22 X2021.01.24 X2021.01.25 X2021.01.26 X2021.01.27
# 1           0          NA           0           1           0           0           0           0
# 2           1          NA           0           1           1           0           0           0
# 3           1          NA           2           1           0           0           0           0
# 4           1          NA           2           0           0           2           0           0
# 5           0          NA           1           0           0           0           0           0
# 6           1          NA           1           0           0           0           0           0
#   Total_Flies
# 1           1
# 2           3
# 3           4
# 4           5
# 5           1
# 6           2

#Convert G2 from wide format into long format. 
G2_long <- G2 %>% pivot_longer(c(X2021.01.15:X2021.01.27), 
                               names_to = "Dates", values_to = "Counts")
#Remove X before dates
G2_long$Dates <- gsub("X","",G2_long$Dates)
#Rplace . with -
G2_long$Dates <- gsub("\\.","-",G2_long$Dates)


#Get the difference in time (Days) between setup and scan date
Diff_Date<- as.numeric(difftime(G2_long$Dates, G2_long$Setup_Date, tz = "UTC"))
G2_long$Diff_Date <- Diff_Date
head(G2_long)
#   Replicate Population Treatment  Bottle Setup_Date Experimentor Embryos Sex   Total_Flies Dates      Counts Diff_Date
#       <int> <chr>      <chr>       <int> <chr>      <chr>          <int> <chr>       <int> <chr>       <int>     <dbl>
# 1         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-15      0        12
# 2         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-16      0        13
# 3         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-18      0        15
# 4         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-19     NA        16
# 5         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-20      0        17
# 6         9 Zinc       25mM_ZnCl2      1 2021-01-03 SJM               25 M               1 2021-01-22      1        19

#Replace NA counts with 0. These were days that no flies emerged from the vials. 
G2_long$Counts <- replace_na(G2_long$Counts, 0)

#Get the time of emergence for each fly
#Code is doing: Take the value in the Diff_Date column, and repeat x times, with
#x being the value in the Counts column.
Devl_Rep <- rep(G2_long$Diff_Date, G2_long$Counts)
G2_long$Replicate <- as.character(G2_long$Replicate)
#Doing same as above, but now with different columns being repeated. 
Rep_Rep <- rep(G2_long$Replicate, G2_long$Counts)
Pop_Rep <- rep(G2_long$Population, G2_long$Counts)
sex_Rep <- rep(G2_long$Sex, G2_long$Counts)
#cbind vectors together
Devl <- cbind(Rep_Rep, Pop_Rep)
Devl <- as.data.frame(Devl)
Devl <- cbind(Devl, sex_Rep)
Devl <- cbind(Devl, Devl_Rep)
colnames(Devl) <- c("Replicate", "Population","Sex", "Emergence_Day")
Devl$Treatment <- rep(G2_long$Treatment, G2_long$Counts)
head(Devl)
#   Replicate Population Sex Emergence_Day  Treatment
# 1         9       Zinc   M            19 25mM_ZnCl2
# 2         9       Zinc   M            15 25mM_ZnCl2
# 3         9       Zinc   M            19 25mM_ZnCl2
# 4         9       Zinc   M            21 25mM_ZnCl2
# 5         9       Zinc   M            15 25mM_ZnCl2
# 6         9       Zinc   M            17 25mM_ZnCl2
dim(Devl)
#[1] 699 5
#Each row represents a fly that emerged during the experiment. Each column
#provides identifying information about the fly such as which replicate and 
#population was it from, what its sex was, the day post setup it emerged,
#and on what treatment it developed. 


#Figure 4A Population Development####
Devl$Treatment <- gsub("0", "O", Devl$Treatment)
Devl$Replicate <- paste("R", Devl$Replicate, sep = "")
Devl$Replicate <- factor(Devl$Replicate, levels = c("R9","R10"))
Devl$Population[Devl$Population == "Control"] <- c("Control Non-Selected")
Devl$Population[Devl$Population == "Zinc"] <- c("Zinc Selected")
Devl$Population <- factor(Devl$Population, levels = c("Control Non-Selected","Zinc Selected"))
Devl$Treatment <- factor(Devl$Treatment, levels = c("H2O","25mM_ZnCl2"))

anova_label <- data.frame(Treatment = "H2O", Emergence_Day = 9)
Devl_anova_rep <- aov(lm(Emergence_Day~Replicate*Population*Treatment*Sex,data = Devl))
anova_devl <- anova(Devl_anova_rep)

Subpop_Devl <- 
  Devl%>% ggplot(aes(x=Treatment, y=Emergence_Day))+
  geom_boxplot(aes(fill=Population),outlier.shape = NA, position = position_dodge(1),
               linewidth=1)+
  geom_point(position = position_jitterdodge(jitter.height = 0.3, dodge.width = 1),
             aes(group=Population, shape=Replicate),
             color="black",fill="white",
             size=2)+
  scale_y_reverse()+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"))+
  geom_text(data = anova_label, aes(x=Treatment, y=Emergence_Day),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova_devl[7,5], digits = 4)),
            nudge_x = .5, size=5)+
  ggtitle("A. X-QTL Subpopulations")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2])))+
  ylab("Development Time (Days)")+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.15,.3))
Subpop_Devl

#Development Analysis####
#A note regarding replicates: As mentioned above, due to the pandemic flies were
#not counted every day and replicate 9 and 10 were setup on different days so
#counting flies on day 10 for replicate 9 would be day 9 for replicate 10.

#Compare each population within the same treatment
#Control H2O Treatment
t.test(Devl %>% filter(Treatment == "H2O" & Population == "Zinc Selected") %>% pull(Emergence_Day),
       Devl %>% filter(Treatment == "H2O" & Population == "Control Non-Selected") %>% pull(Emergence_Day))
# data:  ZW$Emergence_Day and CW$Emergence_Day
# t = -2.8494, df = 348.8, p-value = 0.00464
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.37113197 -0.06801558
# sample estimates:
# mean of x mean of y 
# 11.62835  11.84793 

#Conclusion, our zinc selected population has a faster development time on 
#control media then our non-selected population.
# *** This is not quite the most accurate estimate because we waited until many
# of the flies on control media had emerged. 

t.test(Devl %>% filter(Treatment == "25mM_ZnCl2" & Population == "Zinc Selected") %>% pull(Emergence_Day),
       Devl %>% filter(Treatment == "25mM_ZnCl2" & Population == "Control Non-Selected") %>% pull(Emergence_Day))
# data:  ZZ$Emergence_Day and CZ$Emergence_Day
# t = -3.5483, df = 156.67, p-value = 0.0005119
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.6664796 -0.4746139
# sample estimates:
# mean of x mean of y 
# 16.57143  17.64198 

#Conclusion, our zinc selected population develops faster on zinc media then
#the control non-selected population. 

#Because the zinc selected population had a faster development on control media
#we wanted to know if the difference on zinc media was because the flies were 
#inherently faster at developing or if the zinc treatment had a different effect
#on both population. Population*Treatment Effect
Devl_aov <- aov(lm(Emergence_Day~Replicate*Population*Treatment*Sex,data = Devl))
anova(Devl_aov)
# Response: Emergence_Day
#                                     Df Sum Sq Mean Sq   F value    Pr(>F)    
# Replicate                            1  103.3   103.3   59.9639 3.494e-14 ***
# Population                           1    3.5     3.5    2.0467  0.152999    
# Treatment                            1 4128.9  4128.9 2397.0328 < 2.2e-16 ***
# Sex                                  1    1.7     1.7    0.9671  0.325752    
# Replicate:Population                 1   18.3    18.3   10.6158  0.001177 ** 
# Replicate:Treatment                  1    2.6     2.6    1.4813  0.223992    
# Population:Treatment                 1   23.3    23.3   13.5094  0.000256 ***
# Replicate:Sex                        1    1.3     1.3    0.7702  0.380469    
# Population:Sex                       1    4.9     4.9    2.8180  0.093668 .  
# Treatment:Sex                        1    0.7     0.7    0.4206  0.516833    
# Replicate:Population:Treatment       1    4.3     4.3    2.4754  0.116101    
# Replicate:Population:Sex             1    3.1     3.1    1.8283  0.176775    
# Replicate:Treatment:Sex              1    5.1     5.1    2.9393  0.086900 .  
# Population:Treatment:Sex             1    6.2     6.2    3.6130  0.057750 .  
# Replicate:Population:Treatment:Sex   1   13.7    13.7    7.9511  0.004945 ** 
# Residuals                          683 1176.5     1.7             

#A note regarding replicates: As mentioned above, due to the pandemic flies were
#not counted every day and replicate 9 and 10 were setup on different days so
#counting flies on day 10 for replicate 9 would be day 9 for replicate 10.This can
#be accounting for the significant difference between replicates. 

#There is a significant difference between the effect that zinc treatment had on 
#on both populations. When looking at population development figure and the numbers
#from Devl_stats, we can see that both populations had a developmental delay due
#to zinc. However, the zinc selected population had a shorter delay than the 
#control non-selected population. 

#Sex and Development Time####
#Get the mean and confidence interval for development time for 
#population/treatment/sex groups
Devl_Sex_Stats <- summarySE(Devl, measurevar = "Emergence_Day", 
                            groupvars = c("Population","Sex", "Treatment"), 
                            na.rm = TRUE)
Devl_Sex_Stats
#   Population Sex  Treatment   N Emergence_Day        sd         se        ci
# 1    Control   F 25mM_ZnCl2  44      17.81818 2.3944969 0.36098398 0.7279936
# 2    Control   F        H20 106      11.84906 0.9640436 0.09363619 0.1856633
# 3    Control   M 25mM_ZnCl2  37      17.43243 2.0074935 0.33002989 0.6693316
# 4    Control   M        H20 111      11.84685 1.0108340 0.09594412 0.1901387
# 5       Zinc   F 25mM_ZnCl2  89      16.35955 1.6938873 0.17955170 0.3568213
# 6       Zinc   F        H20 116      11.63793 0.5500443 0.05107033 0.1011605
# 7       Zinc   M 25mM_ZnCl2  51      16.94118 2.5409586 0.35580537 0.7146561
# 8       Zinc   M        H20 145      11.62069 0.6673847 0.05542328 0.1095483

DevlSex <- 
  Devl_Sex_Stats %>% ggplot(aes(x=Population, y=Emergence_Day, color=Sex))+
  geom_errorbar(aes(ymin=Emergence_Day-ci, ymax=Emergence_Day+ci), width=0.3, size=1, position = position_dodge(0.35))+
  geom_point(size=3, position = position_dodge(0.35))+theme_classic()+
  ylab("Development Time (Days)")+
  scale_y_reverse()+
  scale_color_manual(values = c("#CC79A7", "#0072B2"))+
  facet_wrap(~Treatment, strip.position = "bottom")+
  theme(legend.position = c(.8, .8), text = element_text(size = 15))
DevlSex

Devl_aov <- aov(lm(Emergence_Day~Replicate*Population*Treatment*Sex,data = Devl))
anova(Devl_aov)
# Response: Emergence_Day
#                          Df Sum Sq Mean Sq   F value    Pr(>F)
# Population                 1    0.8     0.8    0.4370 0.5087733
# Treatment                  1 4180.7  4180.7 2264.2413 < 2.2e-16 ***
# Sex                        1    1.1     1.1    0.6207 0.4310750
# Population:Treatment       1   25.4    25.4   13.7309 0.0002278 ***
# Population:Sex             1    3.5     3.5    1.9041 0.1680688
# Treatment:Sex              1    1.5     1.5    0.7951 0.3728622
# Population:Treatment:Sex   1    8.4     8.4    4.5628 0.0330242 *
# Residuals                691 1275.9     1.8

#Significant interaction between population*treatment*sex. When looking at DevlSex
#graph made above, this looks to be due to difference in how females respond
#to zinc treatment. In the zinc selected population females tend to experience
#less of a delay than females in the control non-selected population, but the males
#from both populations seem to experience similar delays. 
#This conslusion seems to be supported below when compare each group to each other.
#With the only significant t-test happening when comparing non-selected females on 
#zinc media with zinc selected females on zinc media.

#Comparing development zinc selected males and females on zinc media
#There is not a significant difference
ZZF <- Devl %>% filter(Population == "Zinc Selected" & Treatment == "25mM_ZnCl2" & Sex == "F")
ZZM <- Devl %>% filter(Population == "Zinc Selected" & Treatment == "25mM_ZnCl2" & Sex == "M")
t.test(ZZF$Emergence_Day, ZZM$Emergence_Day)
#p=0.1486

#Comparing development of zinc selected females and non-selected females on zinc media
#Here we do have significant difference
ZCF <- Devl %>% filter(Population == "Control Non-Selected" & Treatment == "25mM_ZnCl2" & Sex == "F")
t.test(ZCF$Emergence_Day, ZZF$Emergence_Day)
#p-value 0.0005814

#Comparing development of non-selected males and females on zinc media
#There is no significant difference
ZCM <- Devl %>% filter(Population == "Control Non-Selected" & Treatment == "25mM_ZnCl2" & Sex == "M")
t.test(ZCM$Emergence_Day, ZCF$Emergence_Day)
#p-value = 0.432

#Comparing development zinc selected males and non0selected males on zinc media
#There is no significant differnece. 
t.test(ZCM$Emergence_Day, ZZM$Emergence_Day)
#p-value 0.3143
