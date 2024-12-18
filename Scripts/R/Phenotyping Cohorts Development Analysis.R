#Introduction####
#Title: Phenotype Population II Development Analysis
#Purpose: To analyze developmental resistance of phenotyping population II, 
#with emergence and development time as phenotypes. Create supplementary material
#figure 3 graphs. 
#Created: 11/18/22
#Last Edited: 7/10/23
#Packages:
library(tidyverse)
library(Rmisc)

#For this analysis you will need "PP-II Development Analysis Scans.txt", and
# "PP-II G2 eggs Barcodes 11012022.txt". 

#Read in Scans####
AllScans <- read.table("PP-II Development Analysis Scans.txt",header = TRUE, 
                       colClasses = "character")
head(AllScans)
#       BC7 Males Females DateDDMMYY TimeHHMM ScanYYMMDDHHMM ScanExperimenter
# 1 2800601     0       0     111122     1137           1137              022
# 2 2800606     0       0     111122     1137           1137              022
# 3 2800611     0       0     111122     1137           1137              022
# 4 2800616     0       0     111122     1137           1137              022
# 5 2800621     0       0     111122     1137           1137              022
# 6 2800626     1       2     111122     1139           1139              022

#This table has all scans from PP-II development assay all together.

#BC7: unique 7 digit barcode for each vial.
#Males and Females: Number of each sex counted that scan.
#Date DDMMYY: Date in terms of day, month, year.
#TimeHHMM: Time in terms of hours and minutes.
#ScanYYMMDDHHM and ScanExperimenter are columns are "left over" columns made
#from the function that reads in the scan files. These have no significance for
#analysis. 

#Read in Barcodes and Merge with All Scans####
Barcode <- read.table("PP-II G2 eggs Barcodes 11012022.txt",header = TRUE)
colnames(Barcode)[1] <- c("BC7")
head(Barcode)
#       BC7          Population   Treatment Number_Eggs Setup_Date
# 1 2800601       Zinc_Selected  25mM_ZnCl2          50 2022-11-01
# 2 2800606       Zinc_Selected         H2O          50 2022-11-01
# 3 2800611       Zinc_Selected   2mM_CuSO4          50 2022-11-01
# 4 2800616       Zinc_Selected 0.2mM_CdCl2          50 2022-11-01
# 5 2800621 Control_Nonselected  25mM_ZnCl2          50 2022-11-01
# 6 2800626 Control_Nonselected         H2O          50 2022-11-01
#Has identifying information for each vial. 
#Population: Eggs are collected from either a zinc selected or control
#non-selected population.
#Treatment: 4 different treatments that eggs could have been put on,
#25mM ZnCl2, 2mM CuSO4, 0.2mM CdCl2 or a control H2O media. 
#Number_Eggs: number of eggs put in the vial.
#Setup_Date: Date that eggs were collected and put in vials. Experiment
#was setup over 4 days. Each population had 10 vials for each treatment created
#each day. 2 population * 4 treatments * 10 vials * 4 setup days = 320 vials

#Convert All_Scans date to format that is in Barcodes yyyy-mm-dd
AllScans$Formatted_Date <- strptime(AllScans$DateDDMMYY, "%d%m%y")

#Combine All_Scans with Barcode
Barcode$BC7 <- as.character(Barcode$BC7)
AllScans$BC7 <- as.character(AllScans$BC7)
AllData <- right_join(Barcode,AllScans,by=c("BC7"))
head(AllData)
#       BC7    Population  Treatment Number_Eggs Setup_Date Males Females DateDDMMYY TimeHHMM ScanYYMMDDHHMM
# 1 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       0     111122     1137           1137
# 2 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       1     201122     1153           1153
# 3 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       0     211122     1102           1102
# 4 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       1     231122     1053           1053
# 5 2800606 Zinc_Selected        H2O          50 2022-11-01     0       0     111122     1137           1137
# 6 2800606 Zinc_Selected        H2O          50 2022-11-01     4       5     131122     1100           1100
#   ScanExperimenter Formatted_Date
# 1               22     2022-11-11
# 2               22     2022-11-20
# 3               22     2022-11-21
# 4               22     2022-11-23
# 5               22     2022-11-11
# 6               22     2022-11-13

#Remove potential F1 flies####
#It is possible, especially in control treatment vials, for flies that emerge
#to mate and lay eggs, and these eggs can then develop on that control media 
#and emerge. 
#Pull out scans where there is at least one male and one female
F1 <- AllData %>% filter(Males >= 1 & Females >= 1)

#Calculate the time between setup and the scan and minus two days. 
#Reason why do minus 2 is because flies were counted in the morning so we assume 
#eggs could have been laid previous day
F1$Time_Between <- as.numeric(difftime(F1$Formatted_Date, F1$Setup_Date, tz = "UTC"))-2
#Add Time_Between value to the scan date to get the date that vial scans should 
#be tossed.
F1$Toss_Date <- as.Date(F1$Formatted_Date) + F1$Time_Between

#If there are duplicate scans that has more then 1 male and female, we want to
#use the earliest toss date.
#Arrange F1 by Formatted_Date so earliest scan is first.
F1 <- F1 %>% arrange(Formatted_Date)
#Keep only the first iteration of a barcode, which will have the earliest toss date
F1 <- F1 %>% distinct(BC7, .keep_all = TRUE)
dim(F1)
#[1] 236 14
#The row number should either be the same number as Barcode or less.

#Isolate only the barcode number and the toss date 
Toss_Date <- F1[c(1,14)]
head(Toss_Date)
#       BC7  Toss_Date
# 1 2800626 2022-11-19
# 2 2800649 2022-11-19
# 3 2800646 2022-11-21
# 4 2800607 2022-11-21
# 5 2800648 2022-11-21
# 6 2800668 2022-11-21

#Join Toss_Date and AllData by barcode
AllData <- full_join(AllData, Toss_Date, by=c("BC7"))

#Not all vials have a toss date, will replace the toss date with 30 days
#post setup
#Calculate what is the date 30 days post setup.
AllData$Day30 <- as.Date(AllData$Setup_Date) + 30
#If there is a NA in toss date column, replace that with day 30 date.
AllData$Toss_Date <- ifelse(is.na(AllData$Toss_Date),as.character(AllData$Day30), 
                            as.character(AllData$Toss_Date))
AllData$Toss_Date <- as.Date(AllData$Toss_Date)

#Remove scans that are on the toss date or after. 
KeepDat <- AllData %>% filter(Formatted_Date < Toss_Date)
head(KeepDat)
#       BC7    Population  Treatment Number_Eggs Setup_Date Males Females DateDDMMYY TimeHHMM
# 1 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       0     111122     1137
# 2 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       1     201122     1153
# 3 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       0     211122     1102
# 4 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       1     231122     1053
# 5 2800606 Zinc_Selected        H2O          50 2022-11-01     0       0     111122     1137
# 6 2800606 Zinc_Selected        H2O          50 2022-11-01     4       5     131122     1100
# ScanYYMMDDHHMM ScanExperimenter Formatted_Date  Toss_Date      Day30
# 1           1137              022     2022-11-11 2022-12-13 2022-12-01
# 2           1153              022     2022-11-20 2022-12-13 2022-12-01
# 3           1102              022     2022-11-21 2022-12-13 2022-12-01
# 4           1053              022     2022-11-23 2022-12-13 2022-12-01
# 5           1137              022     2022-11-11 2022-11-23 2022-12-01
# 6           1100              022     2022-11-13 2022-11-23 2022-12-01
dim(KeepDat)
#[1] 1889 14

RmvDat <- AllData %>% filter(Formatted_Date >= Toss_Date)
dim(RmvDat)
#[1] 47 14
#The row number of KeepDat and RmvDat should equal number of rows in AllData

#Remove unnecessary columns
KeepDat <- KeepDat[-c(8,9,10,11)]
head(KeepDat)
#       BC7    Population  Treatment Number_Eggs Setup_Date Males Females Formatted_Date  Toss_Date      Day30
# 1 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       0     2022-11-11 2022-12-13 2022-12-01
# 2 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     0       1     2022-11-20 2022-12-13 2022-12-01
# 3 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       0     2022-11-21 2022-12-13 2022-12-01
# 4 2800601 Zinc_Selected 25mM_ZnCl2          50 2022-11-01     1       1     2022-11-23 2022-12-13 2022-12-01
# 5 2800606 Zinc_Selected        H2O          50 2022-11-01     0       0     2022-11-11 2022-11-23 2022-12-01
# 6 2800606 Zinc_Selected        H2O          50 2022-11-01     4       5     2022-11-13 2022-11-23 2022-12-01

#Create Emergence Table####
KeepDat$Males <- as.numeric(KeepDat$Males)
KeepDat$Females <- as.numeric(KeepDat$Females)

#Aggregate sex counts by barcodes
Emergence <- aggregate(Males~BC7+Population+Treatment+Number_Eggs+Setup_Date
                     ,KeepDat, sum)
Female_Count <- aggregate(Females~BC7,KeepDat, sum)
Emergence <- merge(Emergence, Female_Count, by=c("BC7"))
Emergence$Total_Flies <- Emergence$Males+Emergence$Females
Emergence$Emergence <- Emergence$Total_Flies/Emergence$Number_Eggs
head(Emergence)
#       BC7          Population  Treatment Number_Eggs Setup_Date Males Females Total_Flies Emergence
# 1 2800301       Zinc_Selected 25mM_ZnCl2          50 2022-11-04     1       5           6      0.12
# 2 2800302 Control_Nonselected 25mM_ZnCl2          50 2022-11-04     0       6           6      0.12
# 3 2800306       Zinc_Selected        H2O          50 2022-11-04    27      16          43      0.86
# 4 2800307 Control_Nonselected        H2O          50 2022-11-04    31      17          48      0.96
# 5 2800311       Zinc_Selected  2mM_CuSO4          50 2022-11-04    16       5          21      0.42
# 6 2800312 Control_Nonselected  2mM_CuSO4          50 2022-11-04    19      22          41      0.82

#Create Development Time Table####
#Calculate difference in days between scan and setup date
Diff_Date<- as.numeric(difftime(KeepDat$Formatted_Date, KeepDat$Setup_Date, tz = "UTC"))
KeepDat$Diff_Date <- Diff_Date
#Create a table where each fly's emergence is a data point
#Create a TOD, Population, Treatment, Replicate, and Generation vector

#Code is doing: Take the value in the TOD column, and repeat x times, with
#x being the value in the Chang_Dead column.
Male_Rep <- rep(KeepDat$Diff_Date, KeepDat$Males)
Male_Rep_BC7 <- rep(KeepDat$BC7, KeepDat$Males)
Male_Rep <- cbind(Male_Rep_BC7, Male_Rep)
Male_Rep <- as.data.frame(Male_Rep)
Male_Rep$Male_Rep <- as.numeric(Male_Rep$Male_Rep)

Female_Rep <- rep(KeepDat$Diff_Date, KeepDat$Females)
Female_Rep_BC7 <- rep(KeepDat$BC7, KeepDat$Females)
Female_Rep <- cbind(Female_Rep_BC7, Female_Rep)
Female_Rep <- as.data.frame(Female_Rep)
Female_Rep$Female_Rep <- as.numeric(Female_Rep$Female_Rep)

Male_Rep$Sex <- c("M")
Female_Rep$Sex <- c("F")

#Combine Female_Rep and Male_Rep
colnames(Female_Rep) <- c("BC7","Diff_Date","Sex")
colnames(Male_Rep) <- c("BC7","Diff_Date","Sex")

All_rep <- rbind(Female_Rep, Male_Rep)
All_rep <- merge(Barcode, All_rep, by=c("BC7"))

head(All_rep)
#       BC7    Population  Treatment Number_Eggs Setup_Date Diff_Date Sex
# 1 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        20   F
# 2 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        17   F
# 3 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        18   F
# 4 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        20   M
# 5 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        20   F
# 6 2800301 Zinc_Selected 25mM_ZnCl2          50 2022-11-04        21   F
dim(All_rep)
# [1] 6933    7

#Clean Up Data####
#2800347 is a Control Population H2O Treatment vial from 4th day of setup that 
#had no flies emerged. After looking at the vial (12/2/22), SJM and KMH suspected that too
#many eggs were put into vial causing the vial to become very wet. It was so wet
#nothing could emerge from the vial. Will remove this vial from our analysis. 
Emergence %>%  filter(Total_Flies == 0)
#        BC7          Population   Treatment Number_Eggs Setup_Date Males Females Total_Flies Emergence
# 1  2800342 Control_Nonselected  25mM_ZnCl2          50 2022-11-04     0       0           0         0
# 2  2800347 Control_Nonselected         H2O          50 2022-11-04     0       0           0         0
# 3  2800357 Control_Nonselected 0.2mM_CdCl2          50 2022-11-04     0       0           0         0
# 4  2800377       Zinc_Selected 0.2mM_CdCl2          50 2022-11-04     0       0           0         0
# 5  2800401 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 6  2800405 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 7  2800416 Control_Nonselected 0.2mM_CdCl2          50 2022-11-03     0       0           0         0
# 8  2800422 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 9  2800425       Zinc_Selected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 10 2800439 Control_Nonselected 0.2mM_CdCl2          50 2022-11-03     0       0           0         0
# 11 2800441 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 12 2800443 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 13 2800445 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 14 2800481 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 15 2800485 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 16 2800500 Control_Nonselected 0.2mM_CdCl2          50 2022-11-03     0       0           0         0
# 17 2800501       Zinc_Selected  25mM_ZnCl2          50 2022-11-02     0       0           0         0
# 18 2800565 Control_Nonselected  25mM_ZnCl2          50 2022-11-03     0       0           0         0
# 19 2800602 Control_Nonselected  25mM_ZnCl2          50 2022-11-01     0       0           0         0
# 20 2800623 Control_Nonselected  25mM_ZnCl2          50 2022-11-01     0       0           0         0
# 21 2800661 Control_Nonselected  25mM_ZnCl2          50 2022-11-01     0       0           0         0
# 22 2800682 Control_Nonselected  25mM_ZnCl2          50 2022-11-01     0       0           0         0
#2800347 is the only control H2O treatment vial not to have any emergence. 

Emergence <- Emergence %>% filter(BC7 != "2800347")
#Do not need to remove from All_rep because no flies emerged from all rep. 

#Check for vials with more then 50 flies emerging
Emergence %>% filter(Total_Flies > 50)
#       BC7          Population Treatment Number_Eggs Setup_Date Males Females Total_Flies Emergence
# 1 2800430       Zinc_Selected       H2O          50 2022-11-03    34      20          54      1.08
# 2 2800568 Control_Nonselected       H2O          50 2022-11-02    25      30          55      1.10
#These two vials had more then 50 eggs put into vial. Will remove from analysis. 

Emergence <- Emergence %>% filter(BC7 != "2800430" & BC7 != "2800568")
All_rep <- All_rep %>% filter(BC7 != "2800430" & BC7 != "2800568")

#Population Emergence Analysis####
#Looking at total population not accounting for sex. 
Population_Treatment <- aggregate(Total_Flies~Population+Treatment, Emergence, sum)
Population_Treatment
#            Population   Treatment Total_Flies
# 1 Control_Nonselected 0.2mM_CdCl2         633
# 2       Zinc_Selected 0.2mM_CdCl2         217
# 3 Control_Nonselected  25mM_ZnCl2         118
# 4       Zinc_Selected  25mM_ZnCl2         236
# 5 Control_Nonselected   2mM_CuSO4        1084
# 6       Zinc_Selected   2mM_CuSO4        1187
# 7 Control_Nonselected         H2O        1652
# 8       Zinc_Selected         H2O        1697

#Emergence using Fisher exact test. 
#Using the counts gives us more power. 
#Testing the control non-selected population and selected population on zinc media
#and H2O media
fisher.test(matrix(c(Population_Treatment[7,3],Population_Treatment[3,3],
                     Population_Treatment[8,3],Population_Treatment[4,3]),ncol = 2))
# p-value = 9.382e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.537432 2.475109
# sample estimates:
# odds ratio 
# 1.94661 

#Testing Control and Copper Treatment
fisher.test(matrix(c(Population_Treatment[7,3],Population_Treatment[5,3],
                     Population_Treatment[8,3],Population_Treatment[6,3]),ncol = 2))
# p-value = 0.2424
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.9567565 1.1876943
# sample estimates:
# odds ratio 
# 1.065952 

#Testing Control and Cadmium Treatment
fisher.test(matrix(c(Population_Treatment[7,3],Population_Treatment[1,3],
                     Population_Treatment[8,3],Population_Treatment[2,3]),ncol = 2))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.2806688 0.3959772
# sample estimates:
# odds ratio 
# 0.3338053 

#When comparing our populations there is a significant difference on zinc where 
#Selected population has a greater emergence then nonselected population, 
#but also on cadmium where the nonselected population has greater emergence then
#zinc selected. 

#Sup Fig 3A & Fig 3C####
Emergence$Treatment <- factor(Emergence$Treatment, levels = c("H2O", "25mM_ZnCl2",
                                                              "0.2mM_CdCl2","2mM_CuSO4"))

Zinc_emerg <- Emergence %>% filter(Treatment == "H2O" | Treatment == "25mM_ZnCl2")
Cohort_Emerg_aov_Zinc <- aov(lm(Emergence~Population*Treatment*Sex,data = Zinc_emerg))
Cohort_anova_Emerg_Zinc <- anova(Cohort_Emerg_aov_Zinc)
# Response: Emergence
#                       Df  Sum Sq Mean Sq   F value    Pr(>F)    
# Population             1  0.0495  0.0495    8.7341  0.003618 ** 
# Treatment              1 23.9418 23.9418 4222.6495 < 2.2e-16 ***
# Population:Treatment   1  0.0332  0.0332    5.8629  0.016633 *  
# Residuals            153  0.8675  0.0057
Cohort_zinc_emerg_label <- data.frame(Treatment = "25mM_ZnCl2", Emergence = 1)

Cu_emerg <- Emergence %>% filter(Treatment == "H2O" | Treatment == "2mM_CuSO4")
Cohort_Emerg_aov_Cu <- aov(lm(Emergence~Population*Treatment,data = Cu_emerg))
Cohort_anova_Emerg_Cu <- anova(Cohort_Emerg_aov_Cu)
# Response: Emergence
#                       Df Sum Sq Mean Sq  F value Pr(>F)    
# Population             1 0.0321  0.0321   1.5919 0.2090    
# Treatment              1 3.5770  3.5770 177.5160 <2e-16 ***
# Population:Treatment   1 0.0252  0.0252   1.2520 0.2649    
# Residuals            153 3.0830  0.0202
Cohort_Cu_emerg_label <- data.frame(Treatment = "2mM_CuSO4", Emergence = 1)

Cd_emerg <- Emergence %>% filter(Treatment == "H2O" | Treatment == "0.2mM_CdCl2")
Cohort_Emerg_aov_Cd <- aov(lm(Emergence~Population*Treatment,data = Cd_emerg))
Cohort_anova_Emerg_Cd <- anova(Cohort_Emerg_aov_Cd)
# Response: Emergence
#                       Df  Sum Sq Mean Sq F value    Pr(>F)    
# Population             1  0.4031  0.4031  19.835 1.620e-05 ***
# Treatment              1 16.9898 16.9898 836.073 < 2.2e-16 ***
# Population:Treatment   1  0.4275  0.4275  21.039 9.304e-06 ***
# Residuals            153  3.1091  0.0203
Cohort_Cd_emerg_label <- data.frame(Treatment = "0.2mM_CdCl2", Emergence = 1)

Cohort_Emergence_all_metals <- Emergence %>% ggplot(aes(x=Treatment, y=Emergence))+
  geom_boxplot(aes(fill=Population), outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             aes(group=Population,shape=Population),
             color="black",fill="white",
             size=3,show.legend = FALSE)+
  scale_shape_manual(values = c(21,21))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"), labels = c("Control Non-Selected","Zinc Selected"))+
  ggtitle("A. Phenotyping Cohorts")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2]),
                            expression("0.2mM CdCl"[2]),expression("2mM CuSO"[4])))+
  geom_text(data = Cohort_zinc_emerg_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ",round(Cohort_anova_Emerg_Zinc[3,5], digits = 3)), size=5)+
  geom_text(data = Cohort_Cu_emerg_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ",round(Cohort_anova_Emerg_Cu[3,5], digits = 2)),size=5)+
  geom_text(data = Cohort_Cd_emerg_label, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ",round(Cohort_anova_Emerg_Cd[3,5], digits = 6)), size=5)+
  theme_classic()+
  theme(text = element_text(size = 15))
Cohort_Emergence_all_metals

Cohort_zinc_emerg_label_only <- data.frame(Treatment = "H2O", Emergence = 1.1)

Cohort_Emergence_zinc <- Emergence %>% filter(Treatment %in% c("H2O", "25mM_ZnCl2")) %>% ggplot(aes(x=Treatment, y=Emergence))+
  geom_boxplot(aes(fill=Population), outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             aes(group=Population,shape=Population),
             color="black",fill="white",
             size=3,show.legend = FALSE)+
  scale_shape_manual(values = c(21,21))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"), labels = c("Control Non-Selected","Zinc Selected"))+
  ggtitle("C. Phenotyping-Only Cohorts")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2])))+
  geom_text(data = Cohort_zinc_emerg_label_only, aes(x=Treatment, y=Emergence),
            label=paste("Pop*Trt ANOVA P-value = ",round(Cohort_anova_Emerg_Zinc[3,5], digits = 4)), size=5,nudge_x = 0.5)+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.15, .2))
Cohort_Emergence_zinc


#Sex and Emergence Analysis####
Pop_Treat_Sex_Male <- aggregate(Males~Population+Treatment, Emergence, sum)
Pop_Treat_Sex_Female <- aggregate(Females~Population+Treatment, Emergence, sum)
Pop_Treat_Sex <- merge(Pop_Treat_Sex_Female, Pop_Treat_Sex_Male, by=c("Population", "Treatment"))
Pop_Treat_Sex
#            Population   Treatment Females Males
# 1 Control_Nonselected 0.2mM_CdCl2     308   325
# 2 Control_Nonselected  25mM_ZnCl2      86    32
# 3 Control_Nonselected   2mM_CuSO4     535   549
# 4 Control_Nonselected         H2O     860   822
# 5       Zinc_Selected 0.2mM_CdCl2     101   116
# 6       Zinc_Selected  25mM_ZnCl2     144    92
# 7       Zinc_Selected   2mM_CuSO4     561   626
# 8       Zinc_Selected         H2O     889   828

#Compare Within Population H2O Treat M:F ratio to other heavy metal M:F ratio 
#Zinc Selected Populations H2O vs zinc treat
fisher.test(matrix(c(Pop_Treat_Sex[8,3], Pop_Treat_Sex[6,3], Pop_Treat_Sex[8,4], 
                     Pop_Treat_Sex[6,4]), ncol = 2))
# p-value = 0.005324
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.5016876 0.8931894
# sample estimates:
# odds ratio 
# 0.6706664 

#Zinc Selected Population H2O vs Cadmium treat
fisher.test(matrix(c(Pop_Treat_Sex[8,3], Pop_Treat_Sex[5,3], Pop_Treat_Sex[8,4], 
                     Pop_Treat_Sex[5,4]), ncol = 2))
# p-value = 0.2202
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.8991669 1.6175903
# sample estimates:
# odds ratio 
# 1.205265 

#Zinc Selected Population H2O vs Copper treatment
fisher.test(matrix(c(Pop_Treat_Sex[8,3], Pop_Treat_Sex[7,3], Pop_Treat_Sex[8,4], 
                     Pop_Treat_Sex[7,4]), ncol = 2))
# p-value = 0.0375
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.006696 1.362425
# sample estimates:
# odds ratio 
# 1.171051
#Boneferoni (0.05/6 = 0.00833333)

#Non-selected Population H2O vs Zinc treatment
fisher.test(matrix(c(Pop_Treat_Sex[4,3], Pop_Treat_Sex[2,3], Pop_Treat_Sex[4,4],
                     Pop_Treat_Sex[2,4]), ncol = 2))
# p-value = 1.464e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.2395055 0.5774771
# sample estimates:
# odds ratio 
# 0.3759093 

#Non-selected Population H2O vs Cadmium Treat
fisher.test(matrix(c(Pop_Treat_Sex[4,3], Pop_Treat_Sex[1,3], Pop_Treat_Sex[4,4], 
                     Pop_Treat_Sex[1,4]),ncol = 2))
# p-value = 0.513
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.8832415 1.2854063
# sample estimates:
# odds ratio 
# 1.06541

#Non-selected Population H2O vs Copper Treat
fisher.test(matrix(c(Pop_Treat_Sex[4,3], Pop_Treat_Sex[3,3], Pop_Treat_Sex[4,4], 
                     Pop_Treat_Sex[3,4]), ncol = 2))
# p-value = 0.6673
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.8862893 1.2113696
# sample estimates:
# odds ratio 
# 1.036159

#Compare within treatment
#Zinc selected and non-selected population on zinc media
fisher.test(matrix(c(Pop_Treat_Sex[2,3], Pop_Treat_Sex[6,3], Pop_Treat_Sex[2,4], 
                     Pop_Treat_Sex[6,4]), ncol = 2))
# p-value = 0.03322
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.035031 2.883385
# sample estimates:
# odds ratio 
# 1.714436 

#Zinc selected and non-selected population on copper treatment
fisher.test(matrix(c(Pop_Treat_Sex[3,3], Pop_Treat_Sex[7,3], Pop_Treat_Sex[3,4], 
                     Pop_Treat_Sex[7,4]), ncol = 2))
# p-value = 0.3337
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.9190056 1.2866893
# sample estimates:
# odds ratio 
# 1.087368 

#Zinc selected and non-selected population on cadmium treatment
fisher.test(matrix(c(Pop_Treat_Sex[1,3], Pop_Treat_Sex[5,3], Pop_Treat_Sex[1,4],
                     Pop_Treat_Sex[5,4]), ncol = 2))
# p-value = 0.6368
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.7898366 1.5011352
# sample estimates:
# odds ratio 
# 1.088326

#Zinc selected and non-selected population on H2O
fisher.test(matrix(c(Pop_Treat_Sex[4,3], Pop_Treat_Sex[8,3], Pop_Treat_Sex[4,4],
                     Pop_Treat_Sex[8,4]), ncol = 2))
# p-value = 0.5804
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.8382001 1.1043040
# sample estimates:
# odds ratio 
# 0.9620822 

#Summary
#Both population are significant (p < 0.05) for zinc treatment where females are
#at a greater frequency on zinc media. 
#Zinc selected population is also significant (p < 0.05) for copper. 
#When comparing different populations on same treatment, there is a significant 
#difference between males and females on zinc media, but not on the other treatments. 

#Development Analysis####
Zinc_Devl <- All_rep %>% filter(Treatment == "H2O" | Treatment == "25mM_ZnCl2")
Zinc_aov_devl <- aov(lm(Diff_Date~Population*Treatment*Sex, data = Zinc_Devl))
anova(Zinc_aov_devl)
# Response: Diff_Date
#                            Df  Sum Sq Mean Sq   F value    Pr(>F)    
# Population                  1     5.2     5.2    1.7683 0.1836690    
# Treatment                   1 16004.8 16004.8 5394.4777 < 2.2e-16 ***
# Sex                         1    20.8    20.8    6.9985 0.0081923 ** 
# Population:Treatment        1     2.7     2.7    0.8987 0.3432032    
# Population:Sex              1     4.6     4.6    1.5580 0.2120361    
# Treatment:Sex               1    36.2    36.2   12.2049 0.0004822 ***
# Population:Treatment:Sex    1     0.1     0.1    0.0230 0.8793922    
# Residuals                3695 10962.6     3.0    

Cd_Devl <- All_rep %>% filter(Treatment == "H2O" | Treatment == "0.2mM_CdCl2")
Cd_aov_devl <- aov(lm(Diff_Date~Population*Treatment*Sex, data = Cd_Devl))
anova(Cd_aov_devl)
#Response: Diff_Date
#Df  Sum Sq Mean Sq   F value  Pr(>F)    
#Population                  1   406.6   406.6  129.7097 < 2e-16 ***
#Treatment                   1  7310.6  7310.6 2331.8878 < 2e-16 ***
#Sex                         1     4.0     4.0    1.2645 0.26087    
#Population:Treatment        1   873.2   873.2  278.5279 < 2e-16 ***
#Population:Sex              1     0.1     0.1    0.0374 0.84665    
#Treatment:Sex               1     4.6     4.6    1.4738 0.22481    
#Population:Treatment:Sex    1    15.2    15.2    4.8605 0.02753 *  
#Residuals                4191 13139.1     3.1

Cu_Devl <- All_rep %>% filter(Treatment == "H2O" | Treatment == "2mM_CuSO4")
Cu_aov_devl <- aov(lm(Diff_Date~Population*Treatment*Sex, data = Cu_Devl))
anova(Cu_aov_devl)
#Response: Diff_Date
#Df  Sum Sq Mean Sq   F value    Pr(>F)    
#Population                  1     6.9     6.9    1.3250 0.2497522    
#Treatment                   1 22531.7 22531.7 4342.8183 < 2.2e-16 ***
#Sex                         1    10.1    10.1    1.9539 0.1622175    
#Population:Treatment        1   262.1   262.1   50.5169 1.329e-12 ***
#Population:Sex              1     0.2     0.2    0.0328 0.8562952    
#Treatment:Sex               1    64.8    64.8   12.4857 0.0004134 ***
#Population:Treatment:Sex    1     2.8     2.8    0.5393 0.4627428    
#Residuals                5612 29116.5     5.2                        

H2O_Devl <- All_rep %>% filter(Treatment == "H2O")
H2O_aov_devl <- aov(lm(Diff_Date~Population*Sex, data = H2O_Devl))
anova(H2O_aov_devl)
#Response: Diff_Date
#Df Sum Sq Mean Sq F value Pr(>F)    
#Population        1  198.7 198.660 75.1474 <2e-16 ***
#Sex               1    6.1   6.071  2.2965 0.1298    
#Population:Sex    1    2.2   2.175  0.8229 0.3644    
#Residuals      3345 8842.9   2.64

#Sup Fig 3A & Fig4C####
All_rep$Treatment <- factor(All_rep$Treatment, levels = c("H2O","25mM_ZnCl2",
                                                          "0.2mM_CdCl2","2mM_CuSO4"))

Cohort_zinc_devl_label <- data.frame(Treatment = "25mM_ZnCl2", Diff_Date = 7)
Cohort_Cu_devl_label <- data.frame(Treatment = "2mM_CuSO4", Diff_Date = 7)
Cohort_Cd_devl_label <- data.frame(Treatment = "0.2mM_CdCl2", Diff_Date = 7)

Cohort_Development_all_metals <- All_rep%>% ggplot(aes(x=Treatment, y=Diff_Date))+
  geom_point(position = position_jitterdodge(jitter.height = 0.3, dodge.width = 1),
             aes(group=Population,shape=Population),
             color="grey40",fill="white",
             size=1,show.legend = FALSE)+
  geom_boxplot(aes(fill=Population,color=Population),
               outlier.shape = NA, position = position_dodge(1),alpha=.5, linewidth=1)+
  scale_y_reverse()+
  scale_shape_manual(values = c(16,16))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"),labels=c("Control Non-Selected","Zinc Selected"))+
  scale_color_manual(values = c("#00AFBB","#FC4E07"), labels=c("Control Non-Selected","Zinc Selected"))+
  geom_text(data = Cohort_zinc_devl_label, aes(x=Treatment, y=Diff_Date),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Zinc_aov_devl)[4,5], digits = 2)), size=5)+
  geom_text(data = Cohort_Cu_devl_label, aes(x=Treatment, y=Diff_Date),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Cu_aov_devl)[4,5], digits = 12)), size=5)+
  geom_text(data = Cohort_Cd_devl_label, aes(x=Treatment, y=Diff_Date),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Cd_aov_devl)[4,5], digits = 60)), size=5)+
  ggtitle("A. Phenotyping Cohorts")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2]),
                            expression("0.2mM CdCl"[2]),expression("2mM CuSO"[4])))+
  ylab("Development Time (Days)")+
  theme_classic()+
  theme(text = element_text(size = 15))
Cohort_Development_all_metals

Cohort_Development_Zinc <- All_rep %>%filter(Treatment %in% c("H2O", "25mM_ZnCl2")) %>% ggplot(aes(x=Treatment, y=Diff_Date))+
  geom_point(position = position_jitterdodge(jitter.height = 0.3, dodge.width = 1),
             aes(group=Population,shape=Population),
             color="grey40",fill="white",
             size=1,show.legend = FALSE)+
  geom_boxplot(aes(fill=Population,color=Population),
               outlier.shape = NA, position = position_dodge(1),alpha=.5, linewidth=1)+
  scale_y_reverse()+
  scale_shape_manual(values = c(16,16))+
  scale_fill_manual(values = c("#00AFBB","#FC4E07"),labels=c("Control Non-Selected","Zinc Selected"))+
  scale_color_manual(values = c("#00AFBB","#FC4E07"), labels=c("Control Non-Selected","Zinc Selected"))+
  geom_text(data = Cohort_zinc_devl_label, aes(x=Treatment, y=Diff_Date),
            label=paste("Pop*Trt ANOVA P-value = ",round(anova(Zinc_aov_devl)[4,5], digits = 4)), size=5, nudge_x = -.5)+
  ggtitle("C. Phenotyping-Only Cohorts")+
  scale_x_discrete(labels=c(expression("H"[2]*"O"),expression("25mM ZnCl"[2]),
                            expression("0.2mM CdCl"[2]),expression("2mM CuSO"[4])))+
  ylab("Development Time (Days)")+
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(.15, .2))
Cohort_Development_Zinc
