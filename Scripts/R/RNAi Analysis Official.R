#Introduction####
#Title: RNAi Analysis
#Purpose: Analyze RNAi KD data measuring both emergence and development time and 
#identify RNAi hits for each phenotype. Will also created supplemental RNAi figures 
#and Figure 5 in paper. 
#Created:6/8/2022
#Last Edited: 23/20/2024
#Packages
library(tidyverse)
library(Rmisc)
library(gridExtra)
library(ggpubr)

#For this analysis you will need "RNAi Scans.txt".

#Columns labeled effect are referring to a lines genotype and if it is 
#a control or RNAi KD. 

#This code can be broken down into three parts. Preparing data for analysis which
#is section Read in Scans - Control Phenotypes Supplementary Material Table 3 
#Data and Figure 5.At the end of Control Phenotypes you have the option to save 
#the two data tables that the rest of the analysis will be using. 
#The next section, beginning at Individual KD Emergence Phenotype, is the actual 
#analysis and identifying RNAi hits. At
#the beginning of Individual KD Emergence Phenotype you can reload in the two
#data tables and proceed from there. 
#The third section is supplemental emergence trend table and graphing KD 
#phenotypes both as seen in supplemental material and in figure 5.

#All RNAi phenotype graphs follow the same format of code. What each section 
#of code does has been broken down in the Control Phenotypes Supplementary
#Material section where the control genotypes' phenotypes are plotted. 

#Read in Scans####
AllScans <- read.table("RNAi Scans.txt", header = TRUE)
head(AllScans)
#   BC7 Father Balancer         Effect     QTL UAS_Chromosome Setup_Date Treatment Embryo Males Females Formatted_Date  Toss_Date      Day30 Rep
# 1   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       0     2022-06-04 2022-06-24 2022-06-24  R1
# 2   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       0     2022-06-05 2022-06-24 2022-06-24  R1
# 3   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       0     2022-06-06 2022-06-24 2022-06-24  R1
# 4   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       0     2022-06-07 2022-06-24 2022-06-24  R1
# 5   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       6     2022-06-08 2022-06-24 2022-06-24  R1
# 6   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0       0     2022-06-09 2022-06-24 2022-06-24  R1

#This table has every scan from all 5 replicates. 
#*R1 was not scanned, instead counts were recorded in an excel table. R1 data
#has been formatted to match scanner's format. 
#Potential F2 flies has been removed from scans. See PP-II Development for 
#how this is done. 
#BC7 unique number identifier for a vial. *R1 is not 7 digits and goes from 1-234.
#Father: UAS line  
#Balancer: Does phenotype have balancer
#Effect: what gene is KD, or what kind of control is it.
#QTL: which QTL is the gene under
#UAS_chromosome: if the UAS gene is on chr2 or 3.
#Setup_Date: Day eggs were collected
#Treatment: What sort of media control H2O or 10mM ZnCl2 was in vial
#Embryo: Number of eggs put in vial.
#Males: Number of males counted in scan
#Females: Number of females counted in scan.
#Formatted_Date: Date of scan
#Toss_Date: Date that scans from vial should be disregarded.
#Day30: 30 days post setup date.
#Rep: Which replicate (R1:R5) are scans from. 

#Create Emergence Table for each Barcode####
#In this table each row will represent one vial and have its percent emergence
#(Total number of flies/Number of Eggs)
#Get total number of flies that emerged at each scan
AllScans$Males <- as.numeric(AllScans$Males)
AllScans$Females <- as.numeric(AllScans$Females)
AllScans$Total_Flies <- AllScans$Males+AllScans$Females

#Sum up all of vial's scans for Males, Females and Total Flies
Emergence <- aggregate(Males~BC7+Father+Balancer+Effect+QTL+UAS_Chromosome+Setup_Date+
                       Treatment+Embryo+Rep,AllScans, sum)
Females<-aggregate(Females~BC7+Effect,AllScans, sum)
Total_Flies <- aggregate(Total_Flies~BC7+Effect,AllScans,sum)

#Merge total males, total females, and total flies into one table
Emergence <- merge(Emergence, Females, by=c("BC7","Effect"))
Emergence <- merge(Emergence,Total_Flies, by=c("BC7","Effect"))
Emergence$Emergence<-Emergence$Total_Flies/Emergence$Embryo
head(Emergence)
#   BC7         Effect Father Balancer     QTL UAS_Chromosome Setup_Date Treatment Embryo Rep Males Females Total_Flies Emergence
# 1   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1     1      10          11      0.22
# 2  10      RanBP3_KD  40948       No       F           Chr2 2022-05-25       H2O     50  R1     9       8          17      0.34
# 3 100       dj-1B_KD  38999       Sb       G           Chr3 2022-05-26       H2O     50  R1     7       8          15      0.30
# 4 101      pHCl-2_KD  26003       No       G           Chr3 2022-05-26      10mM     50  R1     9       8          17      0.34
# 5 102      pHCl-2_KD  26003       No       G           Chr3 2022-05-26       H2O     50  R1    13      11          24      0.48
# 6 103     CG11318_KD  51792       No       G           Chr3 2022-05-26      10mM     50  R1     3      10          13      0.26

#Create Development Table for each Fly####
#In this table each row will represent a fly that emerged and will have the day
#it emerged post setup, and identifying information. 
#Get the difference in days between scan and setup
Diff_Date<- as.numeric(difftime(AllScans$Formatted_Date, AllScans$Setup_Date))
AllScans$Diff_Date <- Diff_Date
#   BC7 Father Balancer         Effect     QTL UAS_Chromosome Setup_Date Treatment Embryo Males
# 1   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
# 2   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
# 3   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
# 4   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
# 5   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
# 6   1  36304       No Control_attP40 Control           Chr2 2022-05-25      10mM     50     0
#   Females Formatted_Date  Toss_Date      Day30 Rep Total_Flies Diff_Date
# 1       0     2022-06-04 2022-06-24 2022-06-24  R1           0        10
# 2       0     2022-06-05 2022-06-24 2022-06-24  R1           0        11
# 3       0     2022-06-06 2022-06-24 2022-06-24  R1           0        12
# 4       0     2022-06-07 2022-06-24 2022-06-24  R1           0        13
# 5       6     2022-06-08 2022-06-24 2022-06-24  R1           6        14
# 6       0     2022-06-09 2022-06-24 2022-06-24  R1           0        15

#Create a table where each fly's death becomes a data point
#Create an Emergence Day, Barcode and Effect vector and combine. 
#Will create seperate tables for males and females and then combine.

#Code is doing: Take the value in the Diff_Date column (or other columns), 
#and repeat x times, with x being the value in the Males or Females column.
Male_Rep <- rep(AllScans$Diff_Date, AllScans$Males)
Male_Rep_BC7 <- rep(AllScans$BC7, AllScans$Males)
Male_Rep_Effect <- rep(AllScans$Effect, AllScans$Males)
#Get additional barcode information from Emergence Table
Barcode_Info <- Emergence[-c(11:14)]
Barcode_Info$BC7 <- as.character(Barcode_Info$BC7)
#Combine Barcode, effect and diff date into one table
Male_Rep <- cbind(Male_Rep_BC7, Male_Rep_Effect, Male_Rep)
Male_Rep <- as.data.frame(Male_Rep)
colnames(Male_Rep) <- c("BC7", "Effect", "Diff_Date")
#Add additional Barcode information
Male_Rep <- right_join(Barcode_Info, Male_Rep, by=c("BC7", "Effect"))
#Create an informational sex column.
Male_Rep$Sex <- c("M")

#Same steps as above except now females
Female_Rep <- rep(AllScans$Diff_Date, AllScans$Females)
Female_Rep_BC7 <- rep(AllScans$BC7, AllScans$Females)
Female_Rep_Effect <- rep(AllScans$Effect, AllScans$Females)
Female_Rep <- cbind(Female_Rep_BC7, Female_Rep_Effect, Female_Rep)
Female_Rep <- as.data.frame(Female_Rep)
colnames(Female_Rep) <- c("BC7", "Effect", "Diff_Date")
Female_Rep <- right_join(Barcode_Info, Female_Rep, by=c("BC7", "Effect"))
Female_Rep$Sex <- c("F")

#Combine Male_Rep and Female_Rep into one table
All_Rep <- rbind(Female_Rep,Male_Rep)
All_Rep$Diff_Date <- as.numeric(All_Rep$Diff_Date)
head(All_Rep)
#   BC7         Effect Father Balancer     QTL UAS_Chromosome Setup_Date Treatment Embryo Rep Diff_Date Sex
# 1   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
# 2   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
# 3   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
# 4   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
# 5   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
# 6   1 Control_attP40  36304       No Control           Chr2 2022-05-25      10mM     50  R1        14   F
tail(All_Rep)
#       BC7             Effect Father Balancer QTL UAS_Chromosome Setup_Date Treatment Embryo Rep Diff_Date Sex
# 15771  96 Control_Fer1HCH_KD  60000       Sb   G           Chr3 2022-05-26       H2O     25  R1        12   M
# 15772  96 Control_Fer1HCH_KD  60000       Sb   G           Chr3 2022-05-26       H2O     25  R1        12   M
# 15773  96         Fer1HCH_KD  60000       Sb   G           Chr3 2022-05-26       H2O     25  R1        12   M
# 15774  96         Fer1HCH_KD  60000       Sb   G           Chr3 2022-05-26       H2O     25  R1        16   M
# 15775  97 Control_Fer2LCH_KD  44067      Cyo   G           Chr2 2022-05-26      10mM     50  R1        12   M
# 15776  97 Control_Fer2LCH_KD  44067      Cyo   G           Chr2 2022-05-26      10mM     50  R1        13   M

#Clean Up Data####
#Filter out Fer2LCH and Fer1HCH vials. These flies did very poorly and were not
#tested past Replicate 2.
dim(Emergence)
#[1] 1178 14
Emergence <- Emergence %>% filter(Father != "44067" & Father != "60000")
dim(Emergence)
#[1] 1125 14
#Filter out flies with Cyo balancer
Emergence <- Emergence %>% filter(Balancer != "Cyo")
dim(Emergence)
#[1] 1079 14
#filter out vials that have less than 40 and greater than 50
Emergence <- Emergence %>% filter(Embryo >=40 & Embryo <=50)
dim(Emergence)
#[1] 1064 14

dim(All_Rep)
#[1] 15624 12
All_Rep <- All_Rep %>% filter(Father != "44067" & Father != "60000")
dim(All_Rep)
#[1] 15536 12
All_Rep <- All_Rep %>% filter(Balancer != "Cyo")
dim(All_Rep)
#[1] 15236 12
All_Rep <- All_Rep %>% filter(Embryo >=40 & Embryo <=50)
dim(All_Rep)
#[1] 15109 12

#Sup Table 5 Data & Sup Fig 3####
#Isolate control genotype lines
Emergence_Control <- Emergence %>% filter(QTL == "Control")
Devl_Control <- All_Rep %>% filter(QTL == "Control")

#Summarize the Emergence phenotype for control genotypes on each treatment across
#Vials.
#Will get the mean, standard deviation, standard error and 95% confidence interval.
Control_Emerg_Stats <- summarySE(data = Emergence_Control, measurevar = "Emergence",
                                 groupvars = c("Effect", "Father", "Treatment"))
Control_Emerg_Stats
#    Effect Father Treatment  N  Emergence         sd          se         ci
# 1  Control_attP2  36303      10mM 29 0.03034483 0.04161754 0.007728183 0.01583047
# 2  Control_attP2  36303       H2O 29 0.45976522 0.14336306 0.026621851 0.05453239
# 3 Control_attP40  36304      10mM 29 0.22482759 0.12255088 0.022757128 0.04661586
# 4 Control_attP40  36304       H2O 29 0.40275862 0.12115210 0.022497381 0.04608380
# 5            GFP  35786      10mM 29 0.14965517 0.07437755 0.013811564 0.02829171
# 6            GFP  35786       H2O 29 0.47446809 0.11410514 0.021188793 0.04340327
# 7     Luciferase  35788      10mM 28 0.24571429 0.08672770 0.016389994 0.03362949
# 8     Luciferase  35788       H2O 29 0.39264368 0.12404997 0.023035502 0.04718609

#Do the same as above but for development time
Control_Devl_Stats <- summarySE(data = Devl_Control, measurevar = "Diff_Date",
                                groupvars = c("Effect", "Father", "Treatment"))
Control_Devl_Stats
#    Effect Father Treatment   N Diff_Date       sd         se         ci
# 1  Control_attP2  36303      10mM  44  13.97727 2.774373 0.41825248 0.84348651
# 2  Control_attP2  36303       H2O 665  12.40752 1.613461 0.06256734 0.12285366
# 3 Control_attP40  36304      10mM 326  13.82515 2.233337 0.12369312 0.24334025
# 4 Control_attP40  36304       H2O 580  11.66552 1.671208 0.06939314 0.13629296
# 5            GFP  35786      10mM 217  12.52535 1.691633 0.11483555 0.22634173
# 6            GFP  35786       H2O 686  11.02624 1.163639 0.04442793 0.08723128
# 7     Luciferase  35788      10mM 344  12.18023 1.746989 0.09419132 0.18526531
# 8     Luciferase  35788       H2O 567  10.75485 1.171605 0.04920277 0.09664232

#The two tables above are where the data for Supplementary Table 3 come from. 

#Graph the means and CI for each genotype and phenotype
#The layout and process to create these graphs and the other RNAi graphs is 
#almost identical (barring gene names etc.). 
#For the control's emergence and development graph we have broken down what each
#section is doing. 

#Before graphing we change H2O treatment to 1 and 10mM Zn to 2.
#This is to help prevent as much empty white space in our graphs.
Control_Emerg_Stats$Treatment[Control_Emerg_Stats$Treatment != "10mM"] <- 1
Control_Emerg_Stats$Treatment[Control_Emerg_Stats$Treatment == "10mM"] <- 2
Control_Emerg_Stats$Treatment <- as.numeric(Control_Emerg_Stats$Treatment)
#Adjusting names to remove _ and or be more descriptive
Control_Emerg_Stats$Effect[Control_Emerg_Stats$Effect == "Control_attP40"] <- c("attP40 Control")
Control_Emerg_Stats$Effect[Control_Emerg_Stats$Effect == "GFP"] <- c("GFP Control")
Control_Emerg_Stats$Effect[Control_Emerg_Stats$Effect == "Luciferase"] <- c("Luciferase Control")
Control_Emerg_Stats$Effect[Control_Emerg_Stats$Effect == "Control_attP2"] <- c("attP2 Control")


#Plot Emergence phenotypes
Control_Emerg_Plot <- 
  Control_Emerg_Stats %>% 
  ggplot(aes(x=Treatment, y=Emergence, color=Effect))+
  geom_errorbar(aes(ymin=Emergence-ci, ymax=Emergence+ci, linetype=NULL), #plot error bars first
                width=.3, linewidth=1, position = position_dodge(.25))+
  geom_line(aes(group=Effect, linetype=Effect), position = position_dodge(.25), #plot line connecting each treatments means
            linewidth=1)+                                                       #This is to help visualize the Genotype*Treatment interaction
  geom_point(size=3, position = position_dodge(.25))+ #plot a point for each genotype's mean
  scale_color_manual(values = c("#2E86C1","gray70","gray70","gray70"))+
  scale_linetype_manual(values = c("dotdash", "dashed", "dotted", "solid"))+
  ggtitle("RNAi Control Emergence")+
  theme_classic()+
  theme(text = element_text(size=20), legend.position = c(.2,.2))+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels=c(expression("H"[2]*"O"),expression("10mM ZnCl"[2])))+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))
Control_Emerg_Plot

#Before graphing we change H2O treatment to 1 and 10mM Zn to 2.
#This is to help prevent as much empty white space in our graphs.
Control_Devl_Stats$Treatment[Control_Devl_Stats$Treatment != "10mM"] <- 1
Control_Devl_Stats$Treatment[Control_Devl_Stats$Treatment == "10mM"] <- 2
Control_Devl_Stats$Treatment <- as.numeric(Control_Devl_Stats$Treatment)
#Adjusting names to remove _ and or be more descriptive
Control_Devl_Stats$Effect[Control_Devl_Stats$Effect == "Control_attP40"] <- c("attP40 Control")
Control_Devl_Stats$Effect[Control_Devl_Stats$Effect == "GFP"] <- c("GFP Control")
Control_Devl_Stats$Effect[Control_Devl_Stats$Effect == "Luciferase"] <- c("Luciferase Control")
Control_Devl_Stats$Effect[Control_Devl_Stats$Effect == "Control_attP2"] <- c("attP2 Control")

Control_Devl_Plot <- Control_Devl_Stats %>% 
  ggplot(aes(x=Treatment, y=Diff_Date, color=Effect))+
  geom_errorbar(aes(ymin=Diff_Date-ci, ymax=Diff_Date+ci, linetype=NULL), 
                width=.3, linewidth=1, position = position_dodge(.25))+ #plot error bars first
  geom_line(aes(group=Effect, linetype=Effect), position = position_dodge(.25), #plot line connecting each treatments means
            linewidth=1)+                                                       #This is to help visualize the Genotype*Treatment interaction
  geom_point(size=3, position = position_dodge(.25))+ #plot a point for each genotype's mean
  scale_color_manual(values = c("#2E86C1","gray70","gray70","gray70"))+
  scale_linetype_manual(values = c("dotdash", "dashed", "dotted", "solid"))+
  scale_y_reverse()+ #All Development time graphs have an inverse y axis. Closer to top, faster it developed
  ylab("Development Time (Days)")+
  ggtitle("RNAi Control Development")+
  theme_classic()+
  theme(text = element_text(size=20), legend.position = "none")+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels=c(expression("H"[2]*"O"),expression("10mM ZnCl"[2])))+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))
Control_Devl_Plot

#Combine Emergence and Development Plots into one figure
Control_Plots <- grid.arrange(Control_Emerg_Plot, Control_Devl_Plot, ncol=2)

#From what we can visually see in the emergence plot is that the attP2 flies
#were extremely susceptible to 10mM ZnCl2, much more then the other controls.

#This is a good place to save these tables and reload them back in if you need
#to step away from the analysis. 
#This is purely optional, in the next section we will load the data in under
#the exact same names. 
write.table(Emergence, "RNAi Emergence Table.txt")
write.table(All_Rep, "Development Time Table.txt")

#Individual KD Emergence Analysis####

#Load in data if returning to the code. 
#If you already have All_Rep and Emergence tables and have gone all way through 
#Control Phenotype Section in your environement then you do not need to do this.
All_Rep <- read.table("Development Time Table.txt", header = TRUE)
Emergence <- read.table("RNAi Emergence Table.txt", header = TRUE)

#Will be comparing each unique KD against three reaming control GFP, Luciferase,
#attP40
#Some genes had two different UAS KDs tested (MTF-1, Xrp1, and Trpm). Will be testing
#each of these individually so need a unique identifier. Combine Effect_Father
#into one column
Emergence$Effect_Father <- paste(Emergence$Effect, Emergence$Father, sep = "_")
count(Emergence$Effect_Father)
#                       x freq
# 1      ATPsynD_KD_33740   40
# 2          ben_KD_28721   40
# 3      CG11318_KD_51792   40
# 4       CG4496_KD_51428   40
# 5   Control_attP2_36303   58
# 6  Control_attP40_36304   58
# 7        dj-1B_KD_38999   39
# 8             GFP_35786   58
# 9       GluRIB_KD_67843   40
# 10       Hsp22_KD_41709   36
# 11     Luciferase_35788   57
# 12       Mekk1_KD_28587   39
# 13        Mnn1_KD_51862   40
# 14       MTF-1_KD_33381   40
# 15       MTF-1_KD_34094   40
# 16       Ndae1_KD_62177   39
# 17        Ndc1_KD_67275   39
# 18    Nup98-96_KD_28562   41
# 19      Octa2R_KD_50678   40
# 20      pHCl-2_KD_26003   40
# 21      RanBP3_KD_40948   40
# 22        Trpm_KD_51713   40
# 23        Trpm_KD_57871   40
# 24        Xrp1_KD_34521   40
# 25        Xrp1_KD_51054   40
#Each KD has about 40 vials of 50 eggs tested, 
#while the controls have 57 or 58 vials. 

#Individual KD analysis with each control.
#Create a list with each unique UAS KD identified 
Emergence_QTL <- Emergence %>% filter(QTL !="Control")
Effect_Name <- unique(Emergence_QTL$Effect_Father)
Effect_Name <- as.character(Effect_Name)

#What these functions are doing:
#For each unique KD (Effect_Name), do an anova analysis with a control genotype
#(Luciferase, GFP, or attP40 Control depending on which function) measuring Emergence
#and looking at the effects of genotype (Effect), treatment, and replicate. Then
#list the anova tables in one document. 
RNAiTest_Luci <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Luciferase")
    output <- lm(Emergence~Effect*Treatment*Rep, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_GFP <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "GFP")
    output <- lm(Emergence~Effect*Treatment*Rep, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_attP40 <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Control_attP40")
    output <- lm(Emergence~Effect*Treatment*Rep, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_attP2 <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Control_attP2")
    output <- lm(Emergence~Effect*Treatment*Rep, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

#Run the functions 
Emerg_Luci <- RNAiTest_Luci(Emergence, Effect_Name)
Emerg_GFP <- RNAiTest_GFP(Emergence, Effect_Name)
Emerg_P40 <- RNAiTest_attP40(Emergence, Effect_Name)
Emerg_P2 <- RNAiTest_attP2(Emergence, Effect_Name)

#Get identifying information (QTL and UAS_Chromosome) for each KD and 
#match that up in each of the ANOVA Table Lists
KD_Details <- Emergence_QTL[c(15, 5, 6)]
KD_Details <- KD_Details %>% distinct(Effect_Father, .keep_all = TRUE)
Emerg_Luci <- right_join(KD_Details,Emerg_Luci,by=c("Effect_Father"))
Emerg_GFP <- right_join(KD_Details, Emerg_GFP,by=c("Effect_Father"))
Emerg_P40 <- right_join(KD_Details, Emerg_P40, by=c("Effect_Father"))
Emerg_P2 <- right_join(KD_Details, Emerg_P2, by=c("Effect_Father"))

head(Emerg_P40)
#     Effect_Father QTL UAS_Chromosome           Factor Df       Sum Sq      Mean Sq      F value       Pr(>F)
# 1 RanBP3_KD_40948   F           Chr2           Effect  1 2.033779e-06 2.033779e-06 2.115923e-04 9.884314e-01
# 2 RanBP3_KD_40948   F           Chr2        Treatment  1 5.739796e-01 5.739796e-01 5.971626e+01 3.165074e-11
# 3 RanBP3_KD_40948   F           Chr2              Rep  4 2.343481e-01 5.858702e-02 6.095335e+00 2.543256e-04
# 4 RanBP3_KD_40948   F           Chr2 Effect:Treatment  1 2.016876e-02 2.016876e-02 2.098338e+00 1.514676e-01
# 5 RanBP3_KD_40948   F           Chr2       Effect:Rep  4 1.157207e-01 2.893019e-02 3.009868e+00 2.302680e-02
# 6 RanBP3_KD_40948   F           Chr2    Treatment:Rep  4 7.215168e-02 1.803792e-02 1.876647e+00 1.228856e-01

#Summarize the Emergence phenotype for each KD on each treatment across
#Vials.
#Will get the mean, standard deviation, standard error and 95% confidence interval.
Emerg_Stats <- summarySE(data = Emergence, measurevar = "Emergence", 
                         groupvars = c("Effect_Father","QTL","UAS_Chromosome","Treatment"))
head(Emerg_Stats)
#      Effect_Father QTL UAS_Chromosome Treatment  N Emergence         sd          se          ci
# 1 ATPsynD_KD_33740   E           Chr3      10mM 20 0.1733878 0.11132876 0.024893867 0.052103462
# 2 ATPsynD_KD_33740   E           Chr3       H2O 20 0.3650000 0.13200877 0.029518059 0.061782007
# 3     ben_KD_28721   A           Chr3      10mM 20 0.0040000 0.01046297 0.002339591 0.004896819
# 4     ben_KD_28721   A           Chr3       H2O 20 0.2980000 0.16807267 0.037582191 0.078660429
# 5 CG11318_KD_51792   G           Chr3      10mM 20 0.2130000 0.07574889 0.016937967 0.035451573
# 6 CG11318_KD_51792   G           Chr3       H2O 20 0.2930000 0.14426219 0.032258006 0.067516782
#

write.csv(Emerg_Stats, "RNAi_Emergence_stats.csv")

#Individual KD Development Analysis####
#Mirros the Emergence Analysis just done above, except now looking at development
#time
All_Rep$Effect_Father <- paste(All_Rep$Effect, All_Rep$Father, sep = "_")

#What these functions are doing:
#For each unique KD (Effect_Name), do an anova analysis with a control genotype
#(Luciferase, GFP, or attP40 Control depending on which function) measuring 
#Development Time(Diff Date) and looking at the effects of genotype (effect), 
#treatment,replicate and sex. Then list the anova tables in one document.
RNAiTest_Devl_Luci <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Luciferase")
    output <- lm(Diff_Date~Effect*Treatment*Rep*Sex, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_Devl_GFP <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "GFP")
    output <- lm(Diff_Date~Effect*Treatment*Rep*Sex, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_Devl_attP40 <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Control_attP40")
    output <- lm(Diff_Date~Effect*Treatment*Rep*Sex, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}

RNAiTest_Devl_attP2 <- function(DATA, EFFECT_NAME){
  Effectlist = list()
  for (i in EFFECT_NAME) {
    #print(i)
    X <- DATA %>% filter(Effect_Father == i | Effect == "Control_attP2")
    output <- lm(Diff_Date~Effect*Treatment*Rep*Sex, data=X)
    #print(output)
    output_2 <- anova(output)
    aovtable <- as.data.frame(output_2)
    aovtable <- aovtable %>% rownames_to_column("Factor")
    aovtable$Effect_Father <- i
    #print(aovtable)
    Effectlist[[i]] <- aovtable
  }
  BigAOV <- do.call(rbind, Effectlist)
}


#Run the functions
Devl_Luci <- RNAiTest_Devl_Luci(All_Rep, Effect_Name)
Devl_GFP <- RNAiTest_Devl_GFP(All_Rep, Effect_Name)
Devl_P40 <- RNAiTest_Devl_attP40(All_Rep, Effect_Name)
Devl_P2 <- RNAiTest_Devl_attP2(All_Rep, Effect_Name)

#Match up identifying information for each KD
Devl_Luci <- right_join(KD_Details, Devl_Luci, by=c("Effect_Father"))
Devl_GFP <- right_join(KD_Details, Devl_GFP, by=c("Effect_Father"))
Devl_P40 <- right_join(KD_Details, Devl_P40, by=c("Effect_Father"))
Devl_P2 <- right_join(KD_Details, Devl_P2, by=c("Effect_Father"))

head(Devl_P40)
#     Effect_Father QTL UAS_Chromosome           Factor Df     Sum Sq    Mean Sq    F value        Pr(>F)
# 1 RanBP3_KD_40948   F           Chr2           Effect  1  175.51198  175.51198  77.143957  4.271061e-18
# 2 RanBP3_KD_40948   F           Chr2        Treatment  1 1346.96378 1346.96378 592.040017 2.014840e-110
# 3 RanBP3_KD_40948   F           Chr2              Rep  4 1571.50883  392.87721 172.683953 1.252425e-121
# 4 RanBP3_KD_40948   F           Chr2              Sex  1  154.39932  154.39932  67.864170  3.785462e-16
# 5 RanBP3_KD_40948   F           Chr2 Effect:Treatment  1   20.73619   20.73619   9.114315  2.579166e-03
# 6 RanBP3_KD_40948   F           Chr2       Effect:Rep  4   72.35074   18.08769   7.950202  2.434155e-06

#Summarize the Development phenotype for each KD on each treatment across
#Vials.
#Will get the mean, standard deviation, standard error and 95% confidence interval.
Devl_Stats <- summarySE(data = All_Rep, measurevar = "Diff_Date",
                        groupvars = c("Effect_Father","QTL","UAS_Chromosome","Treatment"))
head(Devl_Stats)
#      Effect_Father QTL UAS_Chromosome Treatment   N Diff_Date        sd         se        ci
# 1 ATPsynD_KD_33740   E           Chr3      10mM 173  11.31214 1.6829247 0.12795039 0.2525552
# 2 ATPsynD_KD_33740   E           Chr3       H2O 365  10.76712 1.1687383 0.06117456 0.1202999
# 3     ben_KD_28721   A           Chr3      10mM   4  15.00000 0.8164966 0.40824829 1.2992283
# 4     ben_KD_28721   A           Chr3       H2O 298  12.96309 1.8861483 0.10926164 0.2150251
# 5 CG11318_KD_51792   G           Chr3      10mM 213  12.21127 1.6187186 0.11091277 0.2186331
# 6 CG11318_KD_51792   G           Chr3       H2O 293  10.79863 1.5362285 0.08974743 0.1766338

write.csv(Devl_Stats, "RNAi_Development_stats.csv")

#Emergence RNAi Hits####
###Isolate and Merge Effect*Treatment from each ANOVA####
#Labeling Control for ANOVA Table
Emerg_GFP$Control <- c("GFP")
Emerg_Luci$Control <- c("Luciferase")
Emerg_P40$Control <- c("Control_attP40")
Emerg_P2$Control <- c("Control_attP2")

#We are interested in the Effect*Treatment interaction, and not any of other factors
#Isolate the Interaction between Effect and Treatment 
GFP_EffTreat <- Emerg_GFP %>% filter(Factor == "Effect:Treatment")
Luci_EffTreat <- Emerg_Luci %>% filter(Factor == "Effect:Treatment")
P40_EffTreat <- Emerg_P40 %>% filter(Factor == "Effect:Treatment")
P2_EffTreat <- Emerg_P2 %>% filter(Factor == "Effect:Treatment")


#Merge Effect*Treatment interaction 
Emerg_EffTreatment <- merge(GFP_EffTreat, Luci_EffTreat, by = "Effect_Father")
Emerg_EffTreatment <- merge(Emerg_EffTreatment, P40_EffTreat, by = "Effect_Father")
Emerg_EffTreatment <- merge(Emerg_EffTreatment, P2_EffTreat, by= "Effect_Father")
Emerg_EffTreatment <- Emerg_EffTreatment[c(1,2,3,4,5,6,7,8,9,10,15,16,17,18,19,24,25,26,27,28,
                                           33,34,35,36,37)]
colnames(Emerg_EffTreatment) <- c("Effect_Father","QTL", "UAS_Chromosome", 
                                  "Factor", "DF", "Sum_Sq_GFP", "Mean_Sq_GFP", 
                                  "F_Value_GFP", "P_Value_GFP", "GFP_ID",
                                  "Sum_Sq_Luciferase", "Mean_Sq_Luciferase",
                                  "F_Value_Luciferase","P_Value_Luciferase", 
                                  "Luciferase_ID", "Sum_Sq_attP40", 
                                  "Mean_Sq_attP40", "F_Value_attP40", 
                                  "P_Value_attP40", "attP40_ID",
                                  "Sum_Sq_attP2", 
                                  "Mean_Sq_attP2", "F_Value_attP2", 
                                  "P_Value_attP2", "attP2_ID")
head(Emerg_EffTreatment)
#      Effect_Father QTL UAS_Chromosome           Factor DF  Sum_Sq_GFP Mean_Sq_GFP F_Value_GFP  P_Value_GFP GFP_ID Sum_Sq_Luciferase Mean_Sq_Luciferase
# 1 ATPsynD_KD_33740   E           Chr3 Effect:Treatment  1 0.103869692 0.103869692   13.513409 4.271605e-04    GFP       0.012009453        0.012009453
# 2     ben_KD_28721   A           Chr3 Effect:Treatment  1 0.005619109 0.005619109    0.795834 3.750173e-01    GFP       0.126564171        0.126564171
# 3 CG11318_KD_51792   G           Chr3 Effect:Treatment  1 0.354707654 0.354707654   48.411119 8.550694e-10    GFP       0.025687236        0.025687236
# 4  CG4496_KD_51428   B           Chr3 Effect:Treatment  1 0.094421644 0.094421644   10.840451 1.479818e-03    GFP       0.015174604        0.015174604
# 5   dj-1B_KD_38999   G           Chr3 Effect:Treatment  1 0.106858660 0.106858660   12.991488 5.389199e-04    GFP       0.009752678        0.009752678
# 6  GluRIB_KD_67843   D           Chr2 Effect:Treatment  1 0.093661247 0.093661247   12.741312 5.944670e-04    GFP       0.017134606        0.017134606
#   F_Value_Luciferase P_Value_Luciferase Luciferase_ID Sum_Sq_attP40 Mean_Sq_attP40 F_Value_attP40 P_Value_attP40      attP40_ID Sum_Sq_attP2
# 1           1.335583       0.2513878551    Luciferase  0.0011332111   0.0011332111     0.11117832    0.739701276 Control_attP40    0.3334393
# 2          15.172331       0.0002076551    Luciferase  0.0805981951   0.0805981951     8.43878340    0.004776969 Control_attP40    0.1085351
# 3           2.980482       0.0882845407    Luciferase  0.0572641898   0.0572641898     5.82895715    0.018107362 Control_attP40    0.7226008
# 4           1.509084       0.2230204403    Luciferase  0.0022199096   0.0022199096     0.19745365    0.658014553 Control_attP40    0.3167788
# 5           1.022854       0.3149701061    Luciferase  0.0003260739   0.0003260739     0.03043575    0.861951495 Control_attP40    0.3381874
# 6           1.997134       0.1614283587    Luciferase  0.0030010175   0.0030010175     0.30852873    0.580096690 Control_attP40    0.3053425
#   Mean_Sq_attP2 F_Value_attP2 P_Value_attP2      attP2_ID
# 1     0.3334393      46.11396  1.965294e-09 Control_attP2
# 2     0.1085351      16.47226  1.163685e-04 Control_attP2
# 3     0.7226008     105.30283  3.972125e-16 Control_attP2
# 4     0.3167788      38.25498  2.667145e-08 Control_attP2
# 5     0.3381874      43.41835  4.507111e-09 Control_attP2
# 6     0.3053425      44.19242  3.064998e-09 Control_attP2

###Calculate KD's phenotypic Trend relative to each control####
#Create List for each UAS KD
#Rearrange Emerg_Stats so Controls are the last listed.
Emerg_Stats$QTL <- factor(Emerg_Stats$QTL, levels = c("A","B","C","D","E","F","G","Control"))
Emerg_Stats <- Emerg_Stats %>% arrange(QTL)
#Create List of each KD name
Effect_name <- unique(Emerg_Stats$Effect_Father)

#Function is getting the difference between Zn treatment phenotype and 
#the control treatment phenotype (10mM-H2O)
Emergence_Diff <- function(DATA,EFFECT_NAME){
  Emergdifflist =list()
  for (i in EFFECT_NAME) {
    #print(i)
    x <- DATA %>% filter(Effect_Father == i)
    #print(x)
    y <- x[1,6]-x[2,6]
    #print(y)
    z <- as.data.frame(y)
    #print(z)
    colnames(z) <- c("Emergence_Difference")
    z$Effect_Father <- i
    #print(z)
    Emergdifflist[[i]] <- z
  }
  Emerg_Diff <- do.call(rbind,Emergdifflist)
}

Emerg_Diff <- Emergence_Diff(Emerg_Stats,Effect_name)
head(Emerg_Diff)
#                 Emergence_Difference   Effect_Father
# ben_KD_28721              -0.2940000    ben_KD_28721
# CG4496_KD_51428           -0.1969302 CG4496_KD_51428
# Mnn1_KD_51862             -0.1615000   Mnn1_KD_51862
# Ndae1_KD_62177            -0.1514737  Ndae1_KD_62177
# Trpm_KD_51713             -0.0730000   Trpm_KD_51713
# Trpm_KD_57871             -0.0890000   Trpm_KD_57871

#Function is subtracting the control genotype's phenotype from the
#KD's phenotype. 
Emergence_Trend <- function(DATA, EFFECT_NAME_UAS){
  EmergTrendList = list()
  for (i in EFFECT_NAME_UAS) {
    #print(i)
    list <- NULL
    list$Effect_Father <- i
    list <- as.data.frame(list)
    #print(list)
    GFP <- DATA %>% filter(Effect_Father == i | Effect_Father == "GFP_35786")
    GFP_Diff <- GFP[1,1]-GFP[2,1]
    list$GFP_Diff <- GFP_Diff
    #print(list)
    Luci <- DATA %>% filter(Effect_Father == i | Effect_Father == "Luciferase_35788")
    Luci_Diff <- Luci[1,1]-Luci[2,1]
    list$Luciferase_Diff <- Luci_Diff
    #print(list)
    P40 <- DATA %>% filter(Effect_Father == i | Effect_Father == "Control_attP40_36304")
    P40_Diff <- P40[1,1]-P40[2,1]
    list$attP40_Diff <- P40_Diff
    #print(list)
    P2 <- DATA %>% filter(Effect_Father == i | Effect_Father == "Control_attP2_36303")
    P2_Diff <- P2[1,1]-P2[2,1]
    list$attP2_Diff <- P2_Diff
    EmergTrendList[[i]] <- list
  }
  Emerg_Trend <- do.call(rbind,EmergTrendList)
}

Effect_name_UAS <- Emerg_Stats %>% filter(QTL != "Control")
Effect_name_UAS <- unique(Effect_name_UAS$Effect_Father)

Emerg_Trend <- Emergence_Trend(Emerg_Diff, Effect_name_UAS)
head(Emerg_Trend)
#                   Effect_Father   GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# ben_KD_28721       ben_KD_28721 0.03081291    -0.147070608 -0.11606897  0.1354204
# CG4496_KD_51428 CG4496_KD_51428 0.12788268    -0.050000840 -0.01899920  0.2324902
# Mnn1_KD_51862     Mnn1_KD_51862 0.16331291    -0.014570608  0.01643103  0.2679204
# Ndae1_KD_62177   Ndae1_KD_62177 0.17333923    -0.004544292  0.02645735  0.2779467
# Trpm_KD_51713     Trpm_KD_51713 0.25181291     0.073929392  0.10493103  0.3564204
# Trpm_KD_57871     Trpm_KD_57871 0.23581291     0.057929392  0.08893103  0.3404204

#If value is positive means that KD's emergence decreased less then the control's.

###Boneferroni Analysis####
#We have 21 unique UAS KD's, 4 controls and 2 phenotypes so our Boneferroni threshold will be .05/(21*4*2)
Bonferroni <- .05/(21*4*2)

#merge Effect*Treatment and Trend together
Emerg_EffTreatment <- merge(Emerg_EffTreatment, Emerg_Trend, by = "Effect_Father")

#In order to be considered a RNAi hit, candidates have to been Boneferroni sig 
#for one control and then nominally sig for other two
Emerg_Filtered_Genes <- Emerg_EffTreatment %>% filter(P_Value_GFP < Bonferroni | 
                                                        P_Value_Luciferase < Bonferroni | 
                                                        P_Value_attP40 < Bonferroni |
                                                        P_Value_attP2 < Bonferroni) %>% 
  filter(P_Value_GFP < .05 & P_Value_Luciferase < .05 & P_Value_attP40 < .05 & P_Value_attP2 < 0.05)
Emerg_Filtered_Genes
#     Effect_Father QTL UAS_Chromosome           Factor DF Sum_Sq_GFP Mean_Sq_GFP F_Value_GFP  P_Value_GFP GFP_ID Sum_Sq_Luciferase Mean_Sq_Luciferase
# 1  MTF-1_KD_33381   D           Chr3 Effect:Treatment  1  0.4567797   0.4567797    58.90774 3.490866e-11    GFP        0.05703564         0.05703564
# 2 pHCl-2_KD_26003   G           Chr3 Effect:Treatment  1  0.4308513   0.4308513    45.81175 1.969470e-09    GFP        0.04993861         0.04993861
# 3   Xrp1_KD_34521   E           Chr3 Effect:Treatment  1  0.4119021   0.4119021    53.09565 1.986022e-10    GFP        0.04494683         0.04494683
#   F_Value_Luciferase P_Value_Luciferase Luciferase_ID Sum_Sq_attP40 Mean_Sq_attP40 F_Value_attP40 P_Value_attP40      attP40_ID Sum_Sq_attP2
# 1           6.293730         0.01421837    Luciferase    0.09804952     0.09804952       9.554427    0.002765587 Control_attP40    0.8655338
# 2           4.633696         0.03448197    Luciferase    0.08852882     0.08852882       7.405047    0.008017962 Control_attP40    0.8296997
# 3           4.957726         0.02889727    Luciferase    0.08567571     0.08567571       8.345674    0.005002707 Control_attP40    0.8033212
#   Mean_Sq_attP2 F_Value_attP2 P_Value_attP2      attP2_ID  GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# 1     0.8655338     118.56236  2.546541e-17 Control_attP2 0.2778129      0.09992939    0.130931  0.3824204
# 2     0.8296997      92.25823  7.270232e-15 Control_attP2 0.2698129      0.09192939    0.122931  0.3744204
# 3     0.8033212     109.98500  1.472346e-16 Control_attP2 0.2638129      0.08592939    0.116931  0.3684204

Emerg_Trend[Emerg_Filtered_Genes$Effect_Father,]
#                   Effect_Father  GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# MTF-1_KD_33381   MTF-1_KD_33381 0.2778129      0.09992939    0.130931  0.3824204
# pHCl-2_KD_26003 pHCl-2_KD_26003 0.2698129      0.09192939    0.122931  0.3744204
# Xrp1_KD_34521     Xrp1_KD_34521 0.2638129      0.08592939    0.116931  0.3684204

#Looking at each of these hits trends, the trends all are consistent with each other. 
#Each hit's emergence is decreasing less due to zinc then each of the control's emergence. 

#Development RNAi Hits####
###Isolate and Merge Effect*Treatment from each ANOVA####
#Labeling Control in each table
Devl_GFP$Control <- c("GFP")
Devl_Luci$Control <- c("Luciferase")
Devl_P40$Control <- c("Control_attP40")
Devl_P2$Control <- c("Control_attP2")

#We are interested in the Effect*Treatment interaction, and not any of other factors
#Isolate the Interaction between Effect and Treatment 
GFP_Devl_EffTreat <- Devl_GFP %>% filter(Factor == "Effect:Treatment")
Luci_Devl_EffTreat <- Devl_Luci %>% filter(Factor == "Effect:Treatment")
P40_Devl_EffTreat <- Devl_P40 %>% filter(Factor == "Effect:Treatment")
P2_Devl_EffTreat <- Devl_P2 %>% filter(Factor == "Effect:Treatment")


#Merge Effect*Treatment interaction 
Devl_EffTreatment <- merge(GFP_Devl_EffTreat, Luci_Devl_EffTreat, by = "Effect_Father")
Devl_EffTreatment <- merge(Devl_EffTreatment, P40_Devl_EffTreat, by = "Effect_Father")
Devl_EffTreatment <- merge(Devl_EffTreatment, P2_Devl_EffTreat, by= "Effect_Father")
Devl_EffTreatment <- Devl_EffTreatment[c(1,2,3,4,5,6,7,8,9,10,15,16,17,18,19,24,25,26,27,28,33,34,35,36,37)]
colnames(Devl_EffTreatment) <- c("Effect_Father","QTL", "UAS_Chromosome", 
                                 "Factor", "DF", "Sum_Sq_GFP", "Mean_Sq_GFP", 
                                 "F_Value_GFP", "P_Value_GFP", "GFP_ID",
                                 "Sum_Sq_Luciferase", "Mean_Sq_Luciferase",
                                 "F_Value_Luciferase","P_Value_Luciferase", 
                                 "Luciferase_ID", "Sum_Sq_attP40", 
                                 "Mean_Sq_attP40", "F_Value_attP40", 
                                 "P_Value_attP40", "attP40_ID",
                                 "Sum_Sq_attP2", 
                                 "Mean_Sq_attP2", "F_Value_attP2", 
                                 "P_Value_attP2", "attP2_ID")
head(Devl_EffTreatment)
#      Effect_Father QTL UAS_Chromosome           Factor DF Sum_Sq_GFP Mean_Sq_GFP F_Value_GFP  P_Value_GFP GFP_ID Sum_Sq_Luciferase Mean_Sq_Luciferase
# 1 ATPsynD_KD_33740   E           Chr3 Effect:Treatment  1  26.803308   26.803308  23.4924850 1.393148e-06    GFP        15.7972500         15.7972500
# 2     ben_KD_28721   A           Chr3 Effect:Treatment  1   3.483127    3.483127   2.4410710 1.184645e-01    GFP         3.3000254          3.3000254
# 3 CG11318_KD_51792   G           Chr3 Effect:Treatment  1   1.454782    1.454782   0.9809616 3.221360e-01    GFP         0.6267428          0.6267428
# 4  CG4496_KD_51428   B           Chr3 Effect:Treatment  1   4.154373    4.154373   3.1181161 7.764127e-02    GFP         9.7129541          9.7129541
# 5   dj-1B_KD_38999   G           Chr3 Effect:Treatment  1  19.504397   19.504397  16.5433020 5.029212e-05    GFP        19.3632743         19.3632743
# 6  GluRIB_KD_67843   D           Chr2 Effect:Treatment  1  45.087732   45.087732  28.9264037 8.782054e-08    GFP        48.5606214         48.5606214
#   F_Value_Luciferase P_Value_Luciferase Luciferase_ID Sum_Sq_attP40 Mean_Sq_attP40 F_Value_attP40 P_Value_attP40      attP40_ID Sum_Sq_attP2
# 1         14.2473488       1.669388e-04    Luciferase    95.7519754     95.7519754     55.6758617   1.491202e-13 Control_attP40   13.2027207
# 2          2.3783905       1.232924e-01    Luciferase     0.9214989      0.9214989      0.4348086   5.097680e-01 Control_attP40    1.6730489
# 3          0.4325029       5.108733e-01    Luciferase    33.8302633     33.8302633     16.2953419   5.719200e-05 Control_attP40    3.1731103
# 4          7.4718247       6.344844e-03    Luciferase     7.4651002      7.4651002      3.9165448   4.800624e-02 Control_attP40    0.8830172
# 5         16.9033873       4.168779e-05    Luciferase    97.5767391     97.5767391     54.8275405   2.304011e-13 Control_attP40   16.6681908
# 6         31.8288697       2.028497e-08    Luciferase   129.0096139    129.0096139     60.4948247   1.410537e-14 Control_attP40   29.5777690
#   Mean_Sq_attP2 F_Value_attP2 P_Value_attP2      attP2_ID
# 1    13.2027207     9.3513423  0.0022771833 Control_attP2
# 2     1.6730489     0.9197649  0.3377734957 Control_attP2
# 3     3.1731103     1.7445734  0.1868173455 Control_attP2
# 4     0.8830172     0.5413198  0.4620274732 Control_attP2
# 5    16.6681908    11.3518873  0.0007785292 Control_attP2
# 6    29.5777690    15.6157826  0.0000820479 Control_attP2

###Calculate KD's phenotypic Trend relative to each control####

#Function is getting the difference between Zn treatment phenotype and 
#the control treatment phenotype (10mM-H2O)
Development_Diff <- function(DATA,EFFECT_NAME){
  Devldifflist =list()
  for (i in EFFECT_NAME) {
    #print(i)
    x <- DATA %>% filter(Effect_Father == i)
    #print(x)
    y <- x[1,6]-x[2,6]
    #print(y)
    z <- as.data.frame(y)
    #print(z)
    colnames(z) <- c("Development_Difference")
    z$Effect_Father <- i
    #print(z)
    Devldifflist[[i]] <- z
  }
  Devl_Diff <- do.call(rbind,Devldifflist)
}

Devl_Diff <- Development_Diff(Devl_Stats,Effect_name)
head(Devl_Diff)
#                 Development_Difference   Effect_Father
# ben_KD_28721                  2.036913    ben_KD_28721
# CG4496_KD_51428               1.741774 CG4496_KD_51428
# Mnn1_KD_51862                 1.266406   Mnn1_KD_51862
# Ndae1_KD_62177                2.367242  Ndae1_KD_62177
# Trpm_KD_51713                 1.686445   Trpm_KD_51713
# Trpm_KD_57871                 1.324472   Trpm_KD_57871

Development_Trend <- function(DATA, EFFECT_NAME_UAS){
  DevlTrendList = list()
  for (i in EFFECT_NAME_UAS) {
    #print(i)
    list <- NULL
    list$Effect_Father <- i
    list <- as.data.frame(list)
    #print(list)
    GFP <- DATA %>% filter(Effect_Father == i | Effect_Father == "GFP_35786")
    GFP_Diff <- GFP[1,1]-GFP[2,1]
    list$GFP_Diff <- GFP_Diff
    #print(list)
    Luci <- DATA %>% filter(Effect_Father == i | Effect_Father == "Luciferase_35788")
    Luci_Diff <- Luci[1,1]-Luci[2,1]
    list$Luciferase_Diff <- Luci_Diff
    #print(list)
    P40 <- DATA %>% filter(Effect_Father == i | Effect_Father == "Control_attP40_36304")
    P40_Diff <- P40[1,1]-P40[2,1]
    list$attP40_Diff <- P40_Diff
    #print(list)
    P2 <- DATA %>% filter(Effect_Father == i | Effect_Father == "Control_attP2_36303")
    P2_Diff <- P2[1,1]-P2[2,1]
    list$attP2_Diff <- P2_Diff
    DevlTrendList[[i]] <- list
  }
  Devl_Trend <- do.call(rbind,DevlTrendList)
}

Devl_Trend <- Development_Trend(Devl_Diff, Effect_name_UAS)
head(Devl_Trend)
#                   Effect_Father   GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# ben_KD_28721       ben_KD_28721  0.5378062       0.6115303  -0.1227234  0.4671588
# CG4496_KD_51428 CG4496_KD_51428  0.2426671       0.3163912  -0.4178625  0.1720197
# Mnn1_KD_51862     Mnn1_KD_51862 -0.2327010      -0.1589769  -0.8932305 -0.3033483
# Ndae1_KD_62177   Ndae1_KD_62177  0.8681352       0.9418593   0.2076056  0.7974878
# Trpm_KD_51713     Trpm_KD_51713  0.1873384       0.2610625  -0.4731911  0.1166911
# Trpm_KD_57871     Trpm_KD_57871 -0.1746344      -0.1009103  -0.8351639 -0.2452817

#Positive means the KD had greater developmental delay due to zinc then the 
#control, negative means it had less. 

###Boneferroni Analysis####
#merge Effect*Treatment and Trend together
Devl_EffTreatment <- merge(Devl_EffTreatment, Devl_Trend, by = "Effect_Father")

#In order to be considered a RNAi hit, candidates have to been Boneferroni sig 
#for one control and then nominally sig for other two
Devl_Filtered_Genes <- Devl_EffTreatment %>% filter(P_Value_GFP < Bonferroni | 
                                                      P_Value_Luciferase < Bonferroni | 
                                                      P_Value_attP40 < Bonferroni |
                                                      P_Value_attP2 < Bonferroni) %>% 
  filter(P_Value_GFP < .05 & P_Value_Luciferase < .05 & P_Value_attP40 < .05 & P_Value_attP2 < .05)
Devl_Filtered_Genes
#       Effect_Father QTL UAS_Chromosome           Factor DF Sum_Sq_GFP Mean_Sq_GFP F_Value_GFP  P_Value_GFP GFP_ID Sum_Sq_Luciferase Mean_Sq_Luciferase
# 1  ATPsynD_KD_33740   E           Chr3 Effect:Treatment  1  26.803308   26.803308   23.492485 1.393148e-06    GFP         15.797250          15.797250
# 2    dj-1B_KD_38999   G           Chr3 Effect:Treatment  1  19.504397   19.504397   16.543302 5.029212e-05    GFP         19.363274          19.363274
# 3   GluRIB_KD_67843   D           Chr2 Effect:Treatment  1  45.087732   45.087732   28.926404 8.782054e-08    GFP         48.560621          48.560621
# 4    Hsp22_KD_41709   D           Chr2 Effect:Treatment  1  27.815886   27.815886   21.962442 3.052662e-06    GFP         37.959189          37.959189
# 5    Ndae1_KD_62177   B           Chr2 Effect:Treatment  1  50.712889   50.712889   39.627609 4.111413e-10    GFP         54.619418          54.619418
# 6 Nup98-96_KD_28562   F           Chr3 Effect:Treatment  1  35.032335   35.032335   23.875326 1.134448e-06    GFP         32.032483          32.032483
# 7   Octa2R_KD_50678   E           Chr3 Effect:Treatment  1   6.052651    6.052651    5.219628 2.247754e-02    GFP          4.782295           4.782295
# 8   pHCl-2_KD_26003   G           Chr3 Effect:Treatment  1  36.785311   36.785311   17.566689 2.920284e-05    GFP         48.599616          48.599616
# 9     Xrp1_KD_51054   E           Chr2 Effect:Treatment  1  60.493438   60.493438   38.887920 5.803218e-10    GFP         72.074492          72.074492
#   F_Value_Luciferase P_Value_Luciferase Luciferase_ID Sum_Sq_attP40 Mean_Sq_attP40 F_Value_attP40 P_Value_attP40      attP40_ID Sum_Sq_attP2
# 1          14.247349       1.669388e-04    Luciferase      95.75198       95.75198      55.675862   1.491202e-13 Control_attP40    13.202721
# 2          16.903387       4.168779e-05    Luciferase      97.57674       97.57674      54.827540   2.304011e-13 Control_attP40    16.668191
# 3          31.828870       2.028497e-08    Luciferase     129.00961      129.00961      60.494825   1.410537e-14 Control_attP40    29.577769
# 4          30.768760       3.472593e-08    Luciferase     135.19897      135.19897      73.028521   3.309745e-17 Control_attP40    22.885952
# 5          43.806869       5.155681e-11    Luciferase      13.68551       13.68551       7.336561   6.839582e-03 Control_attP40    18.361084
# 6          22.287535       2.558880e-06    Luciferase     116.40820      116.40820      58.390996   3.747580e-14 Control_attP40    22.355311
# 7           4.238476       3.969578e-02    Luciferase      61.95828       61.95828      35.978879   2.516301e-09 Control_attP40     5.991506
# 8          23.540941       1.337331e-06    Luciferase     172.75826      172.75826      66.819267   5.848469e-16 Control_attP40    22.857836
# 9          47.268938       8.995352e-12    Luciferase     202.60572      202.60572      96.896935   3.313727e-22 Control_attP40    33.371917
#   Mean_Sq_attP2 F_Value_attP2 P_Value_attP2      attP2_ID   GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# 1     13.202721      9.351342  2.277183e-03 Control_attP2 -0.9540911     -0.88036703  -1.6146207 -1.0247385
# 2     16.668191     11.351887  7.785292e-04 Control_attP2 -1.0837137     -1.00998960  -1.7442433 -1.1543611
# 3     29.577769     15.615783  8.204790e-05 Control_attP2 -0.8028720     -0.72914789  -1.4634015 -0.8735193
# 4     22.885952     14.658075  1.355634e-04 Control_attP2 -0.6271078     -0.55338374  -1.2876374 -0.6977552
# 5     18.361084     11.641478  6.665797e-04 Control_attP2  0.8681352      0.94185930   0.2076056  0.7974878
# 6     22.355311     12.719555  3.744979e-04 Control_attP2 -0.6027353     -0.52901125  -1.2632649 -0.6733827
# 7      5.991506      4.205831  4.049335e-02 Control_attP2 -0.1646854     -0.09096134  -0.8252150 -0.2353328
# 8     22.857836      9.342137  2.280047e-03 Control_attP2 -0.6613663     -0.58764221  -1.3218959 -0.7320137
# 9     33.371917     17.896943  2.492361e-05 Control_attP2 -0.8881365     -0.81441241  -1.5486661 -0.9587839

Devl_Trend[Devl_Filtered_Genes$Effect_Father,]
#                       Effect_Father   GFP_Diff Luciferase_Diff attP40_Diff attP2_Diff
# ATPsynD_KD_33740   ATPsynD_KD_33740 -0.9540911     -0.88036703  -1.6146207 -1.0247385
# dj-1B_KD_38999       dj-1B_KD_38999 -1.0837137     -1.00998960  -1.7442433 -1.1543611
# GluRIB_KD_67843     GluRIB_KD_67843 -0.8028720     -0.72914789  -1.4634015 -0.8735193
# Hsp22_KD_41709       Hsp22_KD_41709 -0.6271078     -0.55338374  -1.2876374 -0.6977552
# Ndae1_KD_62177       Ndae1_KD_62177  0.8681352      0.94185930   0.2076056  0.7974878
# Nup98-96_KD_28562 Nup98-96_KD_28562 -0.6027353     -0.52901125  -1.2632649 -0.6733827
# Octa2R_KD_50678     Octa2R_KD_50678 -0.1646854     -0.09096134  -0.8252150 -0.2353328
# pHCl-2_KD_26003     pHCl-2_KD_26003 -0.6613663     -0.58764221  -1.3218959 -0.7320137
# Xrp1_KD_51054         Xrp1_KD_51054 -0.8881365     -0.81441241  -1.5486661 -0.9587839

#All trends are consistent within a candidate gene. 

#Isolating Data for RNAi Trend Table for Paper####
KD_Details <- KD_Details[c(1:2)]
KD_Details_Sep <- KD_Details %>% separate_wider_delim(Effect_Father, "_", names=c("Gene", "UAS_Line"), 
                                                  too_many = "merge")
Emerg_Trend_Data <- Emerg_Trend
#Separate Effect_Father column into a gene and UAS_Line column
Emerg_Trend_Data <- Emerg_Trend %>% 
  separate_wider_delim(Effect_Father, "_", names=c("Gene", "UAS_Line"), 
                       too_many = "merge")
Emerg_Trend_Data <- right_join(KD_Details_Sep, Emerg_Trend_Data, by=c("Gene",
                                                                  "UAS_Line"))
#Remove "KD_" from UAS Line
Emerg_Trend_Data$UAS_Line <- gsub("KD_", "", as.character(Emerg_Trend_Data$UAS_Line))

#Separate Effect_Father column into a gene and UAS_Line column
Development_Trend_Data <- Devl_Trend %>% separate_wider_delim(Effect_Father, "_", names=c("Gene", "UAS_Line"), too_many = "merge")
Development_Trend_Data <- right_join(KD_Details_Sep, Development_Trend_Data, by=c("Gene",
                                                                      "UAS_Line"))
#Remove "KD_" from UAS Line
Development_Trend_Data$UAS_Line <- gsub("KD_", "", as.character(Development_Trend_Data$UAS_Line))

#These two tables are the basis of RNAi Trend Table. P-value information
#comes from ANOVA Tables. 

#RNAi Supplemental Graphs####
Emerg_Stats <- Emerg_Stats %>% separate_wider_delim(Effect_Father, "_", names = c("Gene", "UAS_Line"), too_many = "merge",cols_remove = FALSE)
Emerg_Stats$Gene[Emerg_Stats$Gene == "pHCl-1"] <- c("pHCl-2")
Emerg_Stats$Gene <- paste0(Emerg_Stats$Gene, " KD")

Emerg_Stats$Gene_ID <- Emerg_Stats$Gene

Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "attP2_36303"] <- c("attP2 Control Cross")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "attP40_36304"] <- c("attP40 Control Cross")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "35786"] <- c("GFP Control Cross")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "35788"] <- c("Luciferase Control Cross")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_34094"] <- c("MTF-1 KD 34094")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_33381"] <- c("MTF-1 KD 33381")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_51713"] <- c("Trpm KD 51713")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_57871"] <- c("Trpm KD 57871")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_34521"] <- c("Xrp1 KD 34521")
Emerg_Stats$Gene[Emerg_Stats$UAS_Line == "KD_51054"] <- c("Xrp1 KD 51054")

Emerg_Stats$Treatment[Emerg_Stats$Treatment == "H2O"] <- 1
Emerg_Stats$Treatment[Emerg_Stats$Treatment == "10mM"] <- 2
Emerg_Stats$Treatment <- as.numeric(Emerg_Stats$Treatment)
#    Gene      UAS_Line Effect_Father   QTL   UAS_Chromosome Treatment     N Emergence     sd      se      ci Gene_ID  
#   <chr>     <chr>    <chr>           <fct> <chr>              <dbl> <dbl>     <dbl>  <dbl>   <dbl>   <dbl> <chr>    
# 1 ben KD    KD_28721 ben_KD_28721    A     Chr3                   2    20     0.004 0.0105 0.00234 0.00490 ben KD   
# 2 ben KD    KD_28721 ben_KD_28721    A     Chr3                   1    20     0.298 0.168  0.0376  0.0787  ben KD   
# 3 CG4496 KD KD_51428 CG4496_KD_51428 B     Chr3                   2    20     0.179 0.104  0.0231  0.0484  CG4496 KD
# 4 CG4496 KD KD_51428 CG4496_KD_51428 B     Chr3                   1    20     0.376 0.101  0.0225  0.0471  CG4496 KD
# 5 Mnn1 KD   KD_51862 Mnn1_KD_51862   B     Chr2                   2    20     0.125 0.124  0.0277  0.0579  Mnn1 KD  
# 6 Mnn1 KD   KD_51862 Mnn1_KD_51862   B     Chr2                   1    20     0.286 0.0671 0.0150  0.0314  Mnn1 KD  

Devl_Stats <- Devl_Stats %>% separate_wider_delim(Effect_Father, "_", names = c("Gene", "UAS_Line"), too_many = "merge",cols_remove = FALSE)
Devl_Stats$Gene[Devl_Stats$Gene == "pHCl-1"] <- c("pHCl-2")
Devl_Stats$Gene <- paste0(Devl_Stats$Gene, " KD")

Devl_Stats$Gene_ID <- Devl_Stats$Gene

Devl_Stats$UAS_Line[Devl_Stats$UAS_Line == "35786_SJM"] <- c("35786")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "attP2_36303"] <- c("attP2 Control Cross")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "attP40_36304"] <- c("attP40 Control Cross")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "35786"] <- c("GFP Control Cross")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "35788"] <- c("Luciferase Control Cross")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_34094"] <- c("MTF-1 KD 34094")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_33381"] <- c("MTF-1 KD 33381")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_51713"] <- c("Trpm KD 51713")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_57871"] <- c("Trpm KD 57871")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_34521"] <- c("Xrp1 KD 34521")
Devl_Stats$Gene[Devl_Stats$UAS_Line == "KD_51054"] <- c("Xrp1 KD 51054")

Devl_Stats$Treatment[Devl_Stats$Treatment == "H2O"] <- 1
Devl_Stats$Treatment[Devl_Stats$Treatment == "10mM"] <- 2
Devl_Stats$Treatment <- as.numeric(Devl_Stats$Treatment)
#      Gene       UAS_Line Effect_Father    QTL   UAS_Chromosome Treatment     N Diff_Date    sd     se    ci Gene_ID   
#  <chr>      <chr>    <chr>            <chr> <chr>              <dbl> <dbl>     <dbl> <dbl>  <dbl> <dbl> <chr>     
# 1 ATPsynD KD KD_33740 ATPsynD_KD_33740 E     Chr3                   2   173      11.3 1.68  0.128  0.253 ATPsynD KD
# 2 ATPsynD KD KD_33740 ATPsynD_KD_33740 E     Chr3                   1   365      10.8 1.17  0.0612 0.120 ATPsynD KD
# 3 ben KD     KD_28721 ben_KD_28721     A     Chr3                   2     4      15   0.816 0.408  1.30  ben KD    
# 4 ben KD     KD_28721 ben_KD_28721     A     Chr3                   1   298      13.0 1.89  0.109  0.215 ben KD    
# 5 CG11318 KD KD_51792 CG11318_KD_51792 G     Chr3                   2   213      12.2 1.62  0.111  0.219 CG11318 KD
# 6 CG11318 KD KD_51792 CG11318_KD_51792 G     Chr3                   1   293      10.8 1.54  0.0897 0.177 CG11318 KD

KD_list <- Emerg_Stats %>% filter(QTL != "Control") %>% dplyr::count(Gene_ID) %>% filter(n == 2) %>% pull(Gene_ID)
Dup_list <- Emerg_Stats %>% filter(QTL != "Control") %>% dplyr::count(Gene_ID) %>% filter(n == 4) %>% pull(Gene_ID)

Supplement_Graphs <- function(EMERG_STATS, DEVL_STATS, KD_LIST, DUP_LIST){
  pdf("RNAi Supplemental Graphs.pdf", width = 20, height = 10)
  for (i in KD_LIST) {
    print(i)
    Treat <- c(expression("H"[2]*"O"),expression("10mM ZnCl"[2]))
    Plotted_Gene_Emerg <- EMERG_STATS %>% filter(Gene_ID == i | QTL == "Control")
    Plotted_Gene_Emerg$Gene <- factor(Plotted_Gene_Emerg$Gene, levels = c("attP2 Control Cross","attP40 Control Cross", "GFP Control Cross", "Luciferase Control Cross", i))
    Emerg_Plot <- Plotted_Gene_Emerg %>%
      ggplot(aes(x=Treatment, y=Emergence, color=Gene))+
      geom_errorbar(aes(ymin=Emergence-ci, ymax=Emergence+ci, linetype=NULL), width=.1, linewidth=1, position = position_dodge(.1))+
      geom_line(aes(group=Gene, linetype=Gene), position = position_dodge(.1), linewidth=1)+
      geom_point(size=3, position = position_dodge(.1))+
      scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07"))+
      scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
      ggtitle("Emergence on Zinc")+
      theme_classic()+
      theme(text = element_text(size=20), legend.position = c(.2,.1),
            legend.text = element_text(size = 12.5), legend.title = element_text(size=15))+
      scale_x_continuous(name = "Treatment", breaks = c(1.2,1.8), labels = Treat)+
      labs(color = "Genotype", linetype = "Genotype")+
      guides(fill = guide_legend(keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3.5, keyheight = 1),
             colour = guide_legend(keywidth = 3.5, keyheight = 1, override.aes = list(shape=NA)))
    #print(Emerg_Plot)
    Plotted_Gene_Devl <- DEVL_STATS %>% filter(Gene_ID == i | QTL == "Control")
    Plotted_Gene_Devl$Gene <- factor(Plotted_Gene_Devl$Gene, levels = c("attP2 Control Cross","attP40 Control Cross", "GFP Control Cross", "Luciferase Control Cross", i))
    Devl_Plot <- Plotted_Gene_Devl %>%
      ggplot(aes(x=Treatment, y=Diff_Date, color=Gene))+
      geom_errorbar(aes(ymin=Diff_Date-ci, ymax=Diff_Date+ci, linetype=NULL), width=.1, linewidth=1, position = position_dodge(.1))+
      geom_line(aes(group=Gene, linetype=Gene), position = position_dodge(.1), linewidth=1)+
      geom_point(size=3, position = position_dodge(.1))+
      scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07"))+
      scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
      scale_y_reverse()+
      ylab("Development Time (Days)")+
      ggtitle("Development Time")+
      theme_classic()+
      theme(text = element_text(size=20), legend.position = "none")+
      scale_x_continuous(name = "Treatment", breaks = c(1.2,1.8), labels = Treat)+
      labs(color = "Genotype", linetype = "Genotype")+
      guides(fill = guide_legend(keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3.5, keyheight = 1),
             colour = guide_legend(keywidth = 3.5, keyheight = 1, override.aes = list(shape=NA)))
    
    Graph_arrange <- grid.arrange(Emerg_Plot, Devl_Plot, ncol=2, top=text_grob(i, size =20))
    Graph_arrange
    
  }
  for (j in DUP_LIST) {
    print(j)
    Treat <- c(expression("H"[2]*"O"),expression("10mM ZnCl"[2]))
    Plotted_Gene_Emerg <- EMERG_STATS %>% filter(Gene_ID == j | QTL == "Control")
    Plotted_Gene_Emerg$Gene <- factor(Plotted_Gene_Emerg$Gene, levels = c("attP2 Control Cross","attP40 Control Cross", "GFP Control Cross", 
                                                                          "Luciferase Control Cross", Plotted_Gene_Emerg[[1,1]], Plotted_Gene_Emerg[[3,1]]))
    Emerg_Plot_Dup <- 
      Plotted_Gene_Emerg %>%
      ggplot(aes(x=Treatment, y=Emergence, color=Gene))+
      geom_errorbar(aes(ymin=Emergence-ci, ymax=Emergence+ci, linetype=NULL), width=.1, linewidth=1, position = position_dodge(.1))+
      geom_line(aes(group=Gene, linetype=Gene), position = position_dodge(.1), linewidth=1)+
      geom_point(size=3, position = position_dodge(.1))+
      scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07", "blue1"))+
      scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid", "solid"))+
      ggtitle("Emergence on Zinc")+
      theme_classic()+
      theme(text = element_text(size=20), legend.position = c(.2,.11), 
            legend.text = element_text(size = 12.5), legend.title = element_text(size=15))+
      scale_x_continuous(name = "Treatment", breaks = c(1.2,1.8), labels = Treat)+
      labs(color = "Genotype", linetype = "Genotype")+
      guides(fill = guide_legend(keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3.5, keyheight = 1),
             colour = guide_legend(keywidth = 3.5, keyheight = 1, override.aes = list(shape=NA)))
    # test <- grid.arrange(Emerg_Plot_Dup, ncol=1, top=text_grob(j))
    # test
    Plotted_Gene_Devl <- DEVL_STATS %>% filter(Gene_ID == j | QTL == "Control")
    Plotted_Gene_Devl$Gene <- factor(Plotted_Gene_Devl$Gene, levels = c("attP2 Control Cross","attP40 Control Cross", "GFP Control Cross", 
                                                                        "Luciferase Control Cross", Plotted_Gene_Devl[[9,1]], Plotted_Gene_Devl[[11,1]]))
    Devl_Plot_Dup <- 
      Plotted_Gene_Devl %>%
      ggplot(aes(x=Treatment, y=Diff_Date, color=Gene))+
      geom_errorbar(aes(ymin=Diff_Date-ci, ymax=Diff_Date+ci, linetype=NULL), width=.1, linewidth=1, position = position_dodge(.1))+
      geom_line(aes(group=Gene, linetype=Gene), position = position_dodge(.1), linewidth=1)+
      geom_point(size=3, position = position_dodge(.1))+
      scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07", "blue1"))+
      scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid", "solid"))+
      scale_y_reverse()+
      ylab("Development Time (Days)")+
      ggtitle("Development Time")+
      theme_classic()+
      theme(text = element_text(size=20), legend.position = "none")+
      scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
      labs(color = "Genotype", linetype = "Genotype")+
      guides(fill = guide_legend(keywidth = 1, keyheight = 1),
             linetype=guide_legend(keywidth = 3.5, keyheight = 1),
             colour = guide_legend(keywidth = 3.5, keyheight = 1, override.aes = list(shape=NA)))
    Graph_arrange_Dup <- grid.arrange(Emerg_Plot_Dup, Devl_Plot_Dup, ncol=2, top=text_grob(j, size =20))
    Graph_arrange_Dup
  }
  graphics.off()
}

Supplement_Graphs(Emerg_Stats, Devl_Stats, KD_list, Dup_list)

#Figure 5 RNAi KD of pHCl2, MTF-1 and Ndae1####
#Prepare Emergence Data
#Isolate controls and MTF-1, Ndae1 and pHCl-2 data
Paper_Genes_Emerg <- Emerg_Stats%>% filter(Effect_Father == "Luciferase_35788"| 
                                             Effect_Father == "GFP_35786" |
                                             Effect_Father == "Control_attP40_36304"|
                                             Effect_Father == "Control_attP2_36303" |
                                             Effect_Father == "pHCl-2_KD_26003"|
                                             Effect_Father == "Ndae1_KD_62177"|
                                             Effect_Father == "MTF-1_KD_33381" |
                                             Effect_Father == "MTF-1_KD_34094")
#IF running straight through RNAi Supplementary Materials don't need to run
# Paper_Genes_Emerg$Treatment[Paper_Genes_Emerg$Treatment != "10mM"] <- 1
# Paper_Genes_Emerg$Treatment[Paper_Genes_Emerg$Treatment == "10mM"] <- 2
# Paper_Genes_Emerg$Treatment <- as.numeric(Paper_Genes_Emerg$Treatment)
#
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "Control_attP40_36304"] <- c("attP40 Control")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "Control_attP2_36303"] <- c("attP2 Control")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "GFP_35786"] <- c("GFP Control")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "Luciferase_35788"] <- c("Luciferase Control")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "pHCl-2_KD_26003"] <- c("pHCl-2 KD")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "Ndae1_KD_62177"] <- c("Ndae1 KD")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "MTF-1_KD_34094"] <- c("MTF-1 KD 34094")
Paper_Genes_Emerg$Effect_Father[Paper_Genes_Emerg$Effect_Father == "MTF-1_KD_33381"] <- c("MTF-1 KD 33381")

#Prepare Development Data
Paper_Genes_Devl <- Devl_Stats %>% filter(Effect_Father == "Luciferase_35788"| 
                                            Effect_Father == "GFP_35786" |
                                            Effect_Father == "Control_attP40_36304"|
                                            Effect_Father == "Control_attP2_36303" |
                                            Effect_Father == "pHCl-2_KD_26003"|
                                            Effect_Father == "Ndae1_KD_62177"|
                                            Effect_Father == "MTF-1_KD_34094"|
                                            Effect_Father == "MTF-1_KD_33381")

#IF running straight through RNAi Supplementary Materials don't need to run
# Paper_Genes_Devl$Treatment[Paper_Genes_Devl$Treatment == "10mM"] <- 2
# Paper_Genes_Devl$Treatment[Paper_Genes_Devl$Treatment == "H2O"] <- 1
#
Paper_Genes_Devl$Treatment <- as.numeric(Paper_Genes_Devl$Treatment)
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "Control_attP40_36304"] <- c("attP40 Control")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "Control_attP2_36303"] <- c("attP2 Control")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "GFP_35786"] <- c("GFP Control")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "Luciferase_35788"] <- c("Luciferase Control")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "pHCl-2_KD_26003"] <- c("pHCl-2 KD")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "Ndae1_KD_62177"] <- c("Ndae1 KD")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "MTF-1_KD_34094"] <- c("MTF-1 KD 34094")
Paper_Genes_Devl$Effect_Father[Paper_Genes_Devl$Effect_Father == "MTF-1_KD_33381"] <- c("MTF-1 KD 33381")

Treat <- c("H2O", "10mM ZnCl2")

#Make each gene's graph seperately then combine
###pHCl-2####
pHCl2_Emerg<- 
  Paper_Genes_Emerg %>% filter(Effect_Father != "Ndae1 KD" & 
                                 Effect_Father != "MTF-1 KD 34094" &
                                 Effect_Father != "MTF-1 KD 33381") %>% 
  ggplot(aes(x=Treatment, y=Emergence, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
  ggtitle("E. pHCl-2 KD Emergence on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = c(.225,.25))+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))

pHCl2_Emerg

pHCl2_Devl <- 
  Paper_Genes_Devl %>% filter(Effect_Father != "Ndae1 KD" & 
                                Effect_Father != "MTF-1 KD 34094" &
                                Effect_Father != "MTF-1 KD 33381") %>% 
  ggplot(aes(x=Treatment, y=Diff_Date, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","#FC4E07"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
  scale_y_reverse()+
  ylab("Development Time (Days)")+
  ggtitle("F. pHCl-2 KD Development on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = "none")+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))

pHCl2_Devl

###Ndae1####
Ndae1_Emerg <- 
  Paper_Genes_Emerg %>% filter(Effect_Father != "pHCl-2 KD"& 
                                 Effect_Father != "MTF-1 KD 34094" &
                                 Effect_Father != "MTF-1 KD 33381") %>% 
  ggplot(aes(x=Treatment, y=Emergence, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","mediumblue"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
  ggtitle("C. Ndae1 KD Emergence on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = c(.225,.25))+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))
Ndae1_Emerg

Ndae1_Devl <- 
  Paper_Genes_Devl %>% filter(Effect_Father != "pHCl-2 KD"& 
                                Effect_Father != "MTF-1 KD 34094" &
                                Effect_Father != "MTF-1 KD 33381") %>% 
  ggplot(aes(x=Treatment, y=Diff_Date, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","mediumblue"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid"))+
  scale_y_reverse()+
  ylab("Development Time (Days)")+
  ggtitle("D. Ndae1 KD Development on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = "none")+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))
Ndae1_Devl

###MTF-1####
MTF_Emerg_Plot <- 
  Paper_Genes_Emerg %>% filter(Effect_Father != "pHCl-2 KD" &
                                 Effect_Father != "Ndae1 KD") %>% 
  ggplot(aes(x=Treatment, y=Emergence, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","green4", "darkorange"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid", "solid"))+
  ggtitle("A. MTF-1 KD Emergence on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = c(.225,.3))+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3, keyheight = 1),
         colour = guide_legend(keywidth = 3, keyheight = 1, 
                               override.aes = list(shape=NA)))
MTF_Emerg_Plot

MTF_Devl_Plot <- 
  Paper_Genes_Devl %>% filter(Effect_Father != "pHCl-2 KD" &
                                Effect_Father != "Ndae1 KD") %>%   
  ggplot(aes(x=Treatment, y=Diff_Date, color=Effect_Father))+
  geom_line(aes(group=Effect_Father, linetype=Effect_Father), 
            linewidth=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("gray70","gray70","gray70","gray70","green4", "darkorange"))+
  scale_linetype_manual(values = c("dotdash","dashed", "dotted", "solid", "solid", "solid"))+
  scale_y_reverse()+
  ylab("Development Time (Days)")+
  ggtitle("B. MTF-1 KD Development on Zinc")+
  theme_classic()+
  theme(text = element_text(size=15), legend.position = "none")+
  scale_x_continuous(name = "Treatment", n.breaks = 2, labels = Treat)+
  labs(color = "Genotype", linetype = "Genotype")+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1),
         linetype=guide_legend(keywidth = 3.5, keyheight = 1),
         colour = guide_legend(keywidth = 3.5, keyheight = 1, 
                               override.aes = list(shape=NA)))
MTF_Devl_Plot

###Combine####
RNAi_Paper <- grid.arrange(
  MTF_Emerg_Plot, MTF_Devl_Plot,
  Ndae1_Emerg, Ndae1_Devl,
  pHCl2_Emerg, pHCl2_Devl)


ggsave("Fig5_RNAi.png", RNAi_Paper, width = 12, height = 12.5, units = "in")
