#Introduction####
#Title: X-QTL Replicate Collection Data
#Purpose: This script will take you through our analysis of estimating the
#selective pressure of 25mM ZnCl2 and estimating the number of female
#eggs that were tested per replicate. Will also create table 1 used in the paper.
#Created: 4/12/22
#Last Edited: 6/21/23
#Packages Needed:
library(tidyverse)

#You will not need to make any additional directories as long as "XQTL Bottle Summary.txt"
#is in your working directory. 
#As a note, replicate and block are being used interchangeably in this script. 

#Read in data####
Bottle_Counts <- read.table("XQTL Bottle Summary.txt", fill = TRUE, header = TRUE)
dim(Bottle_Counts)
#[1] 234 32
Bottle_Counts[1:6,1:7]
#   Block Setup_date   Treatment Bottle_Number Zinc_Batch Total_Num_Females Total_Num_Males
# 1     1  10/1/2020 Control_H20             1       <NA>               124             106
# 2     1  10/1/2020 Control_H20             2       <NA>               143             111
# 3     1  10/1/2020 Control_H20             3       <NA>               136             119
# 4     1  10/1/2020 Control_H20             4       <NA>               120             110
# 5     1  10/1/2020   Zinc_25mM             1          A                33              15
# 6     1  10/1/2020   Zinc_25mM             2          A                31              17

#Columns 1:5 give information about each bottle that flies were collected from.
#Block: Designate which replicate bottles are from (1-12)
#Setup_date: Date that eggs were pippetted into bottles
#Treatment: Bottles with instant media and water (Control_H2O) or Bottles with instant media 
#and 25mM ZnCl2 solution (Zinc_25mM)
#Bottle Number: Unique bottle number for bottles within a replicate,treatment group
#Zinc_Batch: 9 batches (A:I) of 25mM ZnCl2 were made through the course of the experiment.
#Control bottles will have an NA.
#Total_Num_Females: Total number of females that emerged and were collected from bottles
#Tot_Num_Males: Total number of males that emerged and were collected from bottles

#Columns 8:32
#Collection Day Information. A day's collection data is divided into 4 columns.

Bottle_Counts[1:6, 8:11]
#   Collection1_Date Collection1_DayAfterSetup Collection1_Num_Females Collection1_Num_Males
# 1       10/12/2020                        11                      30                    19
# 2       10/12/2020                        11                      15                     7
# 3       10/12/2020                        11                      28                    15
# 4       10/12/2020                        11                      11                    12
# 5       10/16/2020                        15                      10                     5
# 6       10/16/2020                        15                      10                     6

#Collection#_Date: Date that flies were collected
#Collection#_DayAfterSetup: Number of days that had passed from setup to collection.
#Collection#_Num_Females: Number of females collected that day.
#Collection#_Num_Males: Number of males collected that day. 
#There were up to 6 collection days for bottles. Bottles that did not have a collection day
#past a certain number, there are NAs in those places. 

#Note: Replicate 4's 4th control bottle has -99 in the Collection2_Num_Males column
#This is because males were tossed before their count was recorded.
#Because our analysis is focused on females, this does not effect our analysis. 

#Group by Replicate and Treatment####
#We will be adding the the counts among the bottles that are in the same replicate
#and treatment. 

#Replace NAs with 0 (cannot use aggregate with NAs)
Bottle_Counts_0 <- Bottle_Counts
Bottle_Counts_0[is.na(Bottle_Counts_0)] <- 0

#Check that zinc batch had no effect within replicate 1.
#Rep1 had two different batches of 25mM ZnCl2 used and want to be sure there is 
#no effect.
Rep1 <- Bottle_Counts_0 %>% filter(Treatment == "Zinc_25mM" & Block == 1)
zinc_batch <- aov(lm(Total_Num_Females~Zinc_Batch, data = Rep1))
anova(zinc_batch)
# Response: Total_Num_Females
#            Df  Sum Sq Mean Sq F value Pr(>F)
# Zinc_Batch  1    4.41   4.408  0.0516 0.8228
# Residuals  18 1536.79  85.377 
#No significant effect so will combine all of replicate 1 zinc vials together

#Aggregate Bottle Counts (total numbers and each collection day) 
#in Bottle_Counts_0 dataframe by Block, Setup Date and Treatment columns
Counts <- aggregate(cbind(Total_Num_Females, Total_Num_Males,
                          Collection1_Num_Females, Collection1_Num_Males, 
                          Collection2_Num_Females, Collection2_Num_Males, 
                          Collection3_Num_Females, Collection3_Num_Males, 
                          Collection4_Num_Females, Collection4_Num_aMales, 
                          Collection5_Num_Females, Collection5_Num_Males,
                          Collect6_Num_Females, Collection6_Num_Males)
                    ~Block+Setup_date+Treatment, 
                    Bottle_Counts_0,sum)
dim(Counts)
#[1] 24 17
Counts[1:3, 1:7]
#   Block Setup_date   Treatment Total_Num_Females Total_Num_Males Collection1_Num_Females Collection1_Num_Males
# 1     1  10/1/2020 Control_H20               523             446                      84                    53
# 2     3 10/19/2020 Control_H20               401             375                     156                   104
# 3     2  10/2/2020 Control_H20               498             497                     110                   115
Counts[15:17, 1:7]
#    Block Setup_date Treatment Total_Num_Females Total_Num_Males Collection1_Num_Females Collection1_Num_Males
# 15     2  10/2/2020 Zinc_25mM               329             192                      56                    23
# 16     4 10/20/2020 Zinc_25mM               260             146                       5                     0
# 17     5 10/24/2020 Zinc_25mM               192              84                      84                    16

#Each row in the counts table has the sum of flies collected within each replicate and treatment group. 

#Estimating female survival#### 
#We can only estimate what percentage of female eggs developed into adults because 
#we do not have exact number of eggs put in each bottle. 
#To estimate this we will calculate the number of adult females that emerged per 
#ul of eggs in the control bottles, and treat this number as 1 to calculate the 
#fraction of eggs on zinc that emerged.

#Create Table with the number of bottles tested, and total number of females 
#collected per replicate/treatment
#Getting number of bottles per replicate/treatment
Bottles_Num <- Bottle_Counts %>% dplyr::count(Block, Treatment)
#Getting the total number of females per replicate/treatment
Total_Females <- Counts[c(1,3,4)]
#Merging the two data frames together
Total_Females$Group <- paste(Total_Females$Treatment, Total_Females$Block, sep = "-")
Bottles_Num$Group <- paste(Bottles_Num$Treatment, Bottles_Num$Block, sep = "-")
#Due to merge() only merging by one group, 
#need combine Replicate and Treatment information into one column
Bottles_Num <- Bottles_Num[c(3,4)]
Total_Females <- Total_Females[c(3,4)]
#Remove the redundant Replicate and Treatment columns
Block_Emergence <- merge(Bottles_Num, Total_Females, by = c("Group"), all = TRUE)
Block_Emergence <- Block_Emergence %>% separate(Group, c("Treatment","Block"), sep = "-", convert = TRUE)
#Merge Bottle_Num and Total_Females and separate Group column into Treatment and 
#Block columns.
dim(Block_Emergence)
#[1] 24 4
head(Block_Emergence)
#     Treatment Block n Total_Num_Females
# 1 Control_H20     1 4               523
# 2 Control_H20    10 4               725
# 3 Control_H20    11 4               638
# 4 Control_H20    12 4               739
# 5 Control_H20     2 4               498
# 6 Control_H20     3 4               401

#Calculating Number of Emerged flies/ul
#Add ul per bottle column
#Replicate 8's control bottles had 48 ul of eggs instead of 24, and zinc bottles had 60ul eggs intead 48
Block_Emergence$ul_per_bottle[!Block_Emergence$Block == 8 & Block_Emergence$Treatment == "Control_H20"] <- 24
Block_Emergence$ul_per_bottle[!Block_Emergence$Block == 8 & Block_Emergence$Treatment == "Zinc_25mM"] <- 48
Block_Emergence$ul_per_bottle[Block_Emergence$Block == 8 & Block_Emergence$Treatment == "Control_H20"] <- 48
Block_Emergence$ul_per_bottle[Block_Emergence$Block == 8 & Block_Emergence$Treatment == "Zinc_25mM"] <- 60
#Calculate total number of ul per replicate/treatment
Block_Emergence$ul_per_block <- Block_Emergence$n*Block_Emergence$ul_per_bottle
#Calculate Number of Emerged females per ul
Block_Emergence$Female_per_ul <- Block_Emergence$Total_Num_Females/Block_Emergence$ul_per_block
head(Block_Emergence)
#     Treatment Block n Total_Num_Females ul_per_bottle ul_per_block Female_per_ul
# 1 Control_H20     1 4               523            24           96      5.447917
# 2 Control_H20    10 4               725            24           96      7.552083
# 3 Control_H20    11 4               638            24           96      6.645833
# 4 Control_H20    12 4               739            24           96      7.697917
# 5 Control_H20     2 4               498            24           96      5.187500
# 6 Control_H20     3 4               401            24           96      4.177083

#Create table that has both control and zinc selected females/ul on same row
Selection_Emergence <- Block_Emergence[c(13:24),]
Selection_Emergence$Control_per_ul <- Block_Emergence[c(1:12),(7)]
head(Selection_Emergence)
#    Treatment Block  n Total_Num_Females ul_per_bottle ul_per_block Female_per_ul Control_per_ul
# 13 Zinc_25mM     1 20               376            48          960     0.3916667       5.447917
# 14 Zinc_25mM    10 12               290            48          576     0.5034722       7.552083
# 15 Zinc_25mM    11 12               311            48          576     0.5399306       6.645833
# 16 Zinc_25mM    12 14               379            48          672     0.5639881       7.697917
# 17 Zinc_25mM     2 17               329            48          816     0.4031863       5.187500
# 18 Zinc_25mM     3 18               299            48          864     0.3460648       4.177083

#Estimate percent of female eggs on zinc media that emerged as adults.
#Divide Zinc females/ul by Control females/ul
Selection_Emergence$PercentEmergence <- Selection_Emergence$Female_per_ul/Selection_Emergence$Control_per_ul*100
#Gives us each replicates estimated percent emergence
# head(Selection_Emergence)
#    Treatment Block  n Total_Num_Females ul_per_bottle ul_per_block Female_per_ul Control_per_ul PercentEmergence
# 13 Zinc_25mM     1 20               376            48          960     0.3916667       5.447917         7.189293
# 14 Zinc_25mM    10 12               290            48          576     0.5034722       7.552083         6.666667
# 15 Zinc_25mM    11 12               311            48          576     0.5399306       6.645833         8.124347
# 16 Zinc_25mM    12 14               379            48          672     0.5639881       7.697917         7.326503
# 17 Zinc_25mM     2 17               329            48          816     0.4031863       5.187500         7.772266
# 18 Zinc_25mM     3 18               299            48          864     0.3460648       4.177083         8.284843
mean(Selection_Emergence$PercentEmergence)
#[1] 7.042155
#Estimated percent emergence on zinc media averaged across all twelve replicates
#This is number we are referring to when we say approximately 7% of females eggs 
#made it to adulthood. 

#Estimating how many females were tested on zinc####
#Using number females/ul on zinc, and the total number of ul of eggs pippetted 
#per block on zinc media, we can estimate how many female eggs were exposed to 
#25mM ZnCl2.
Selection_Emergence$Estimated_Females_Tested <- Selection_Emergence$ul_per_block*Selection_Emergence$Control_per_ul
head(Selection_Emergence)
#    Treatment Block  n Total_Num_Females ul_per_bottle ul_per_block Female_per_ul Control_per_ul PercentEmergence Estimated_Females_Tested
# 13 Zinc_25mM     1 20               376            48          960     0.3916667       5.447917         7.189293                     5230
# 14 Zinc_25mM    10 12               290            48          576     0.5034722       7.552083         6.666667                     4350
# 15 Zinc_25mM    11 12               311            48          576     0.5399306       6.645833         8.124347                     3828
# 16 Zinc_25mM    12 14               379            48          672     0.5639881       7.697917         7.326503                     5173
# 17 Zinc_25mM     2 17               329            48          816     0.4031863       5.187500         7.772266                     4233
# 18 Zinc_25mM     3 18               299            48          864     0.3460648       4.177083         8.284843                     3609
mean(Selection_Emergence$Estimated_Females_Tested)
#[1] 4313.688
#This is number we are reffering to when we say approximately 4000 females eggs 
#were exposed to zinc per replicate. 

#Supplementary Table 1####
#Table presented in the paper is a reorganization/combination of Block_Emergence 
#and Selection_Emergence. Numbers come directly from those two tables. 

#Create table with each block's details on number bottles, ul, emerged females
#Separate Block Emergence into a control treatment and zinc selective treatment
Block_Emergence_C <- Block_Emergence %>% filter(Treatment == "Control_H20")
Block_Emergence_Z <- Block_Emergence %>% filter(Treatment == "Zinc_25mM")
#Change column names to reflect which treatment the numbers are for
colnames(Block_Emergence_C) <- c("Treatment", "Replicate", "Number_Control_Bottles", 
                                 "Total_Control_Females","Control_ul/Bottle",
                                     "Control_ul/Replicate", "Control_Females/ul")
colnames(Block_Emergence_Z) <- c("Treatment", "Replicate", "Number_Zinc_Bottles", 
                                 "Total_Zinc_Females","Zinc_ul/Bottle",
                                  "Zinc_ul/Replicate", "Zinc_Females/ul")
#Merge control and zinc table by block
Replicate_Table <- merge(Block_Emergence_C, Block_Emergence_Z, by = "Replicate")
#Remove treatment columns
Replicate_Table <- Replicate_Table[-c(2,8)]

#Pull out calculated estimates of females tested and selection. 
Selection_Emergence_Estimates <- Selection_Emergence[c(2,10,9)]
colnames(Selection_Emergence_Estimates) <- c("Replicate","Estimated_Females_Tested_Zinc", 
                                             "Estimated_Percent_Females_Emerged_Zinc")
#Merge with Block_Table by Replicate
Replicate_Table <- merge(Replicate_Table, Selection_Emergence_Estimates, by = "Replicate")
Replicate_Table
#    Replicate Number_Control_Bottles Total_Control_Females Control_ul/Bottle Control_ul/Replicate Control_Females/ul Number_Zinc_Bottles
# 1          1                      4                   523                24                   96           5.447917                  20
# 2          2                      4                   498                24                   96           5.187500                  17
# 3          3                      4                   401                24                   96           4.177083                  18
# 4          4                      4                   479                24                   96           4.989583                  21
# 5          5                      4                   437                24                   96           4.552083                  12
# 6          6                      4                   486                24                   96           5.062500                  20
# 7          7                      4                   600                24                   96           6.250000                  20
# 8          8                      4                  1017                48                  192           5.296875                  12
# 9          9                      4                   754                24                   96           7.854167                   8
# 10        10                      4                   725                24                   96           7.552083                  12
# 11        11                      4                   638                24                   96           6.645833                  12
# 12        12                      4                   739                24                   96           7.697917                  14
#    Total_Zinc_Females Zinc_ul/Bottle Zinc_ul/Replicate Zinc_Females/ul Estimated_Females_Tested_Zinc Estimated_Percent_Females_Emerged_Zinc
# 1                 376             48               960       0.3916667                       5230.00                               7.189293
# 2                 329             48               816       0.4031863                       4233.00                               7.772266
# 3                 299             48               864       0.3460648                       3609.00                               8.284843
# 4                 260             48              1008       0.2579365                       5029.50                               5.169500
# 5                 192             48               576       0.3333333                       2622.00                               7.322654
# 6                 190             48               960       0.1979167                       4860.00                               3.909465
# 7                 476             48               960       0.4958333                       6000.00                               7.933333
# 8                 351             60               720       0.4875000                       3813.75                               9.203540
# 9                 169             48               384       0.4401042                       3016.00                               5.603448
# 10                290             48               576       0.5034722                       4350.00                               6.666667
# 11                311             48               576       0.5399306                       3828.00                               8.124347
# 12                379             48               672       0.5639881                       5173.00                               7.326503

Replicate_Table_Rearranged <- Replicate_Table[c(1:2,4,3,6:7,9,8,11:13)]
write.csv(Replicate_Table, file = paste("Sup Table 1 XQTL Block Estimates.csv"), row.names = FALSE)
#Column names were slightly adjusted after being exported out of R. 
