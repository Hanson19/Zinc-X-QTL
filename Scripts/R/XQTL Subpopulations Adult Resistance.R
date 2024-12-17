#Introduction####
#Title: X-QTL Subpopulations Adult Resistance Analysis
#Purpose: To analyze time of death of adult zinc naive female flies derived
#from zinc selected and non-selected flies 
#Created: 6/26/23
#Last Edited: 6/28/23
#Packages:
library(tidyverse)
library(Rmisc)

#For this anlaysis you will need "AllScanData_G1.txt", "AllScanData_G3.txt", 
#"Adult_G1_FC_011221_KMH.DAT", and "Adult_G3_FC_020521_KMH.DAT". 

#Read in Data####
#Read in scans that record number of flies in each vial. 
#There are scans from Generation 1 flies and Generation 3 flies. 
#We will be combining them all into one AllScan file. 

AllScan_G1 <- read.table("AllScanData_G1.txt", header = TRUE, colClasses = "character")
AllScan_G3 <- read.table("AllScanData_G3.txt", header = TRUE, colClasses = "character")
AllScan_G3$SetUpBlockNumeric <- "2"
SubScans <- rbind(AllScan_G1, AllScan_G3)
colnames(SubScans)[10] <- "Pop_Trt"
head(SubScans)
#       BC7 NumDead DateDDMMYY TimeHHMM ScanYYMMDDHHMM ScanExperimenter TimeInHours SetUpExperimenter          SetUpBlock      Pop_Trt SetUpBlockNumeric
# 1 2632751       0     040121     0854     2101040854              KMH           0               KMH F1_Setup_010421_KMH B9_SelZ_TrtC                 1
# 2 2632746       0     040121     0854     2101040854              KMH           0               KMH F1_Setup_010421_KMH B9_SelZ_TrtC                 1
# 3 2632741       0     040121     0854     2101040854              KMH           0               KMH F1_Setup_010421_KMH B9_SelC_TrtC                 1
# 4 2632736       0     040121     0855     2101040855              KMH           0               KMH F1_Setup_010421_KMH B9_SelC_TrtC                 1
# 5 2632797       1     040121     0855     2101040855              KMH           0               KMH F1_Setup_010421_KMH B9_SelC_TrtZ                 1
# 6 2632786       0     040121     0858     2101040858              KMH           0               KMH F1_Setup_010421_KMH B9_SelC_TrtZ                 1

#Subscans has every vial scan that counted the number of dead.
#Each unique BC7 represents 1 vial.
#NumDead is the number of dead flies in vial at time of scanning
#DateDDMMYY is the day the vials were setup
#TimeHHMM is the time that vials were first scanned during setup. 
#ScanYYMMDDHHMM is the date and time of the scan to count the dead was done.
#2101040854 = January 4th, 2021 8:54am
#ScanExperimenter is the person who scanned and counted number of dead
#TimeinHours is the number of hours between Setup (DateDDMMYY TimeHHMM) and
#scan (ScanYYMMDDHHMM)
#SetupExperimentor the person who setup the vials (collecting,putting flies
#in vials, and initial scan)
#SetUpBlock includes information about generation tested (F1 = G1, F3 = G3),
#date of setup, and experimenter who set it up.
#Pop_Trt includes information about which population flies came from (B9_SelZ =
#Replicate 9 Zinc Selected population) and if these adults were on control media 
#(TrtC) or zinc media (TrtZ)
#SetUpBlockNumeric number for experiment block. G1 flies are 1, G3 are 2. 

#Check for Inconsistencies in Scans####
Barcodes <- unique(SubScans$BC7)

# IDENTIFY inconsistent scans
# -> Where a once dead fly appears to come back
#    to life due to experimenter error/confusion
# -> Function adds new column to each element
#    of 'ScansBy7BC'
findinconsistency.func <- function(SCANS, BARCODE){
  z <- NULL
  for (i in BARCODE) {
    DAT <- SCANS %>% filter(BC7 == i)
    
    if(nrow(DAT)==1) {
      InconsistentScan <- 0
    }
    else {
      # EXTRACT 'NumDead' column
      # CREATE new vector
      # --> start with diffs between current #dead and next #dead
      # --> if this is -ve number Y=1, otherwise Y=0
      X <- as.numeric(DAT[,2])
      Y <- 1*(c(diff(X),0)<0)
      
      # LOOP through datapoints (skipping
      # 1st and last)
      for(i in 2:(length(X)-1)) {
        # IF current value of Y is not negative
        if(Y[i]==0) {
          # DO NOTHING
        }
        else {
          # IF Y[i] = -ve AND current & previous
          # values of 'X' are not equal
          if(X[i]!=X[(i-1)]) {
            # DO NOTHING
          }
          # IF not, move the -ve number indicating
          # "1" in Y to _next_ position
          else {
            Y[i] <- 0
            Y[(i+1)] <- 1
          }
        }
      }
      # RETURN Y vector indicating which rows
      # of input 'DAT' should be eliminated
      # APPENDED to original DAT matrix
      InconsistentScan <- Y
    }
    w <- cbind(DAT,InconsistentScan)
   z <- rbind(w,z)
  }
  print(z)
}

SubScans <- findinconsistency.func(SubScans,Barcodes)

FixOneInconsistentScan <- function(XXX) {
  # EXTRACT 'NumDead' and 'InconsistentScan'
  # columns
  cND <- XXX[,"NumDead"]
  cIS <- XXX[,"InconsistentScan"]
  # LOOP through latter
  for(i in 1:length(cIS)) {
    # SET current value of 'cIS'
    cur.cIS <- cIS[i]
    # IF==0, do nothing
    if(cur.cIS==0) {
      # DO NOTHING
    }
    # IF==1, current value of number dead
    # to the value at the next scan 
    else {
      cND[i] <- cND[i+1]
    }
  }
  # ADD the new number dead vector to the
  # original matrix and return
  Output <- cbind(XXX,cND,stringsAsFactors = F)
  colnames(Output)[ncol(Output)] <- "NumDeadFinal"
  return(Output)
}

Gold <- NULL
Valid <- NULL
Manual <- NULL

# LOOP through elements

for(i in Barcodes) {
  # EXTRACT column indicating whether there are
  # inconsistent scans in a single 'ScansBy7BC' element
  ScansBy7BC <- SubScans %>% filter(BC7 == i)
  InconColumn <- as.numeric(ScansBy7BC[,"InconsistentScan"])
  
  # IF zero problems
  # ADD new column to element (identical to 'NumDead')
  # MOVE element to Gold table
  if(sum(InconColumn)==0) {
    ScansBy7BC$NumDeadFinal <- ScansBy7BC$NumDead
    Gold <- rbind(Gold, ScansBy7BC)
  }
  else {
    # IF there is exactly 1 problem, 'solve' it by
    # assuming the _next_ scan is correct
    # MOVE edited element to the Valid table
    if(sum(InconColumn)==1) {
      ScansBy7BC.Validated <- FixOneInconsistentScan(ScansBy7BC)
      Valid <- rbind(Valid, ScansBy7BC.Validated)
    }
    # IF there is >1 problem, move the un-edited
    # element to 'ScansBy7BC.ForManualCorrection'
    # and we will need to solve this manually
    else {
      ScansBy7BC$FixMe <- ScansBy7BC$NumDead
      Manual <- rbind(Manual, ScansBy7BC)
    }
  }
}

dim(SubScans)
#867 12
dim(Gold)
#859 13
dim(Valid)
#8 13
dim(Manual)
#NULL
#Number of rows in Gold, Valid and Manual should equal number of rows in SubScans                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

#Do have scans in Valid which is where there was 1 inconsistency
unique(Valid$BC7)
#"2632395"
#Only one vial had an inconsistency, 2632395

Valid
#          BC7 NumDead DateDDMMYY TimeHHMM ScanYYMMDDHHMM ScanExperimenter      TimeInHours SetUpExperimenter
# 3100 2632395       0     280121     1234     2101281234              KMH 1.23333333333333               KMH
# 901  2632395       0     290121     1111     2101291111              KMH            23.85               KMH
# 1781 2632395       0     300121     1106     2101301106              KMH 47.7666666666667               KMH
# 2661 2632395       1     310121     1112     2101311112              KMH 71.8666666666667               KMH
# 354  2632395       0     010221     1150     2102011150              KMH             96.5               KMH
# 441  2632395       0     020221     1107     2102021107              KMH 119.783333333333               KMH
# 512  2632395       0     030221     1111     2102031111              KMH           143.85               KMH
# 542  2632395       0     040221     1109     2102041109              KMH 167.816666666667               KMH
#               SetUpBlock      Pop_Trt SetUpBlockNumeric InconsistentScan NumDeadFinal
# 3100 F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 901  F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 1781 F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 2661 F3_Setup_012821_KMH B9_SelZ_TrtC                 2                1            0
# 354  F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 441  F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 512  F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0
# 542  F3_Setup_012821_KMH B9_SelZ_TrtC                 2                0            0

#On scan made on January 31st 2021 (4th row) counted the number of dead as 1, but 
#all other scan have 0 counted. The code has corrected that 1 to a 0 in the
#NumDeadFinal column

#rbind Gold and Valid together
corrected_scans <- rbind(Gold, Valid)

#Correct for Initial Dead####
#Occasionally a fly may already be dead during setup and not dead during treatment. 
#We count them on setup day and the days after in scans counting number of death.
#We want to exclude these flies before doing our analysis.

#Put TimeInHours in terms of days
corrected_scans$TimeInHours <- as.numeric(corrected_scans$TimeInHours)
corrected_scans$Day <- round(corrected_scans$TimeInHours/24)

Z <- NULL
for (i in Barcodes) {
  filt_corrected <- corrected_scans %>% filter(BC7 == i)
  Day_0 <- filt_corrected %>% filter(Day == 0)
  filt_corrected$Corrected_Dead <- filt_corrected$NumDeadFinal
  #There is one vial "2632305" that was not scanned on setup.
  #This will cause Day_0 to have no rows and will exclude these vials if we do
  #not have below if statement.
  #We will assume there were no initial dead flies and use the counts scanned
  if (nrow(Day_0) == 0){
    Z <- rbind(filt_corrected, Z)
  }
  #If it does have rows and the value in Day 0's NumDeadFinal column is 0
  #then just use the values that are in the NumDeadFinal column and put them 
  #in new column called Corrected_Dead
  else if (Day_0$NumDeadFinal == 0){
    filt_corrected$Corrected_Dead <- filt_corrected$NumDeadFinal
    Z <- rbind(filt_corrected, Z)
  }
  #If Day 0's NumDeadFinal column is not 0, designated that number as initial dead
  #Substract initial dead from the NumDeadFinal column and have these values go
  #into new column called Corrected_Dead
  else {
    Initial_Dead <- as.numeric(Day_0$NumDeadFinal)
    filt_corrected$NumDeadFinal <- as.numeric(filt_corrected$NumDeadFinal)
    filt_corrected$Corrected_Dead <- filt_corrected$NumDeadFinal-Initial_Dead
    Z <- rbind(filt_corrected,Z)
  }
}
Dead_Correct <- Z
dim(Dead_Correct)
#867 15
dim(corrected_scans)
#867 14
#Dead_Correct should have the same number of rows as corrected_scans. If it has 
#less then lost some scans in loop. 

#Can check to make sure the loops look by checking "2632797" which had a fly dead
#day 0
Dead_Correct %>% filter(BC7 == "2632797")
#         BC7 NumDead DateDDMMYY TimeHHMM ScanYYMMDDHHMM ScanExperimenter TimeInHours SetUpExperimenter
# 5   2632797       1     040121     0855     2101040855              KMH     0.00000               KMH
# 29  2632797       1     050121     1204     2101051204              KMH    27.15000               KMH
# 81  2632797       2     060121     1043     2101061043              KMH    49.80000               KMH
# 129 2632797      11     070121     1050     2101071050              KMH    73.91667               KMH
# 177 2632797      18     080121     1040     2101081040              KMH    97.75000               KMH
# 225 2632797      19     090121     1040     2101091040              KMH   121.75000               KMH
# 272 2632797      20     100121     1053     2101101053              KMH   145.96667               KMH
#              SetUpBlock      Pop_Trt SetUpBlockNumeric InconsistentScan NumDeadFinal Day Corrected_Dead
# 5   F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0            1   0              0
# 29  F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0            1   1              0
# 81  F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0            2   2              1
# 129 F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0           11   3             10
# 177 F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0           18   4             17
# 225 F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0           19   5             18
# 272 F1_Setup_010421_KMH B9_SelC_TrtZ                 1                0           20   6             19

#Inconsistency Final Count####
#After all flies are dead a final scan is made to count total number of flies in a vial. 
#(This is so we can get numbers for control treatment vials)
#Occasionally the last scan of zinc treatment vials and the total flies in vial scan do not match up.
#If this happens we give the last death scan the value of the total flies in vial scan

#Separate Pop_Trt to individual columns so can filter by treatment
Dead_Correct <- Dead_Correct %>% separate(Pop_Trt, c("Replicate", "Population", "Treatment"))

#Read in Total Flies in vial scan
Total_Flies_G3 <- read.csv("Adult_G3_FC_020521_KMH.DAT", header = FALSE)
Total_Flies_G1 <- read.csv("Adult_G1_FC_011221_KMH.DAT", header = FALSE)
Total_Flies <- rbind(Total_Flies_G1, Total_Flies_G3)
colnames(Total_Flies) <- c("BC7", "Total_Dead", "Date", "Time")
Total_Flies$Total_Dead <- as.numeric(Total_Flies$Total_Dead)

#Correct Total Flies in Vial scan by initial dead
Dead_initial <- corrected_scans %>% filter(Day == 0)
Dead_initial <- Dead_initial[,c(1,13)]
colnames(Dead_initial)[2] <- c("Initial_Dead")
Dead_initial$Initial_Dead <- as.numeric(Dead_initial$Initial_Dead)
DontLoose_2632305 <- Total_Flies %>% filter(BC7 == "2632305")
DontLoose_2632305$Initial_Dead <- 0
Total_Flies <- merge(Total_Flies, Dead_initial)
Total_Flies$Total_Flies <- Total_Flies$Total_Dead-Total_Flies$Initial_Dead
DontLoose_2632305$Total_Flies <- DontLoose_2632305$Total_Dead-DontLoose_2632305$Initial_Dead
Total_Flies <- rbind(Total_Flies, DontLoose_2632305)

#Only want to do this for zinc treatment vials
TrtZ_BC <- Dead_Correct %>% filter(Treatment == "TrtZ")
TrtZ_BC <- unique(TrtZ_BC$BC7)
TrtC <- Dead_Correct %>% filter(Treatment == "TrtC")

Final_Count_Fix <- NULL
for (i in TrtZ_BC) {
  x <- Dead_Correct %>% filter(BC7 == i)
  x <- arrange(x, desc(TimeInHours))
  TF <- Total_Flies %>% filter(BC7 == i)
  LC <- x[1,]
  if(LC$Corrected_Dead == TF$Total_Flies){
    Final_Count_Fix <- rbind(Final_Count_Fix, x)
  }
  else{
    x[1,17] <- TF$Total_Flies
    Final_Count_Fix <- rbind(Final_Count_Fix, x)
  }
}
Final_Count_Fix <- rbind(Final_Count_Fix, TrtC)

dim(Final_Count_Fix)
#867 17
dim(Dead_Correct)
#867 17
#Final_Count_Fix and Dead_Correct should match in dimensions

#Change Dead each Day####
#For our analysis we want the time of death (TOD) for each fly. 
#Our scans have the total number of dead flies per day, but we want the change
#in dead per day
Final_Count_Fix$Corrected_Dead <- as.numeric(Final_Count_Fix$Corrected_Dead)
Final_Count_Fix <- Final_Count_Fix %>% arrange(Day)
Dead_Table <- Final_Count_Fix %>% pivot_wider(id_cols = c("BC7", "Replicate", "Population", "Treatment", "SetUpBlockNumeric"),
                                           names_from = "Day", values_from = "Corrected_Dead")
head(Dead_Table)
#   BC7     Replicate Population Treatment SetUpBlockNumeric   `0`   `1`   `2`   `3`   `4`   `5`   `6`   `7`
#   <chr>   <chr>     <chr>      <chr>     <chr>             <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 2632339 B9        SelC       TrtZ      2                     0     0     6    12    18    20    NA    NA
# 2 2632302 B10       SelZ       TrtZ      2                     0     1     2     6    15    18    19    20
# 3 2632385 B9        SelZ       TrtZ      2                     0     0     4    16    19    20    NA    NA
# 4 2632397 B10       SelC       TrtZ      2                     0     0     3    13    19    20    NA    NA
# 5 2632344 B9        SelC       TrtZ      2                     0     0     2    13    19    20    NA    NA
# 6 2632396 B10       SelC       TrtZ      2                     0     0     0     9    18    20    NA    NA

Dead_Table <- merge(Dead_Table,Total_Flies[-c(2:5)], by=c("BC7"))

#Replace NAs
#A vial may not be scanned a day because it was accidentally skipped, but flies
#were still alive. (See "2632326). If that happens we will take the average difference
#between day before missing scan and day after scan as number of flies that died. 
#In this data 2632326 is the only vial this happened to. 

#Control Treatment that Day 7 is NA
CT_Day7 <- Dead_Table %>% filter(BC7 == "2632731" | BC7 == "2632726" | BC7 == "2632716" | BC7 == "2632721")
CT_Day7$`7` <- 0
Dead_Table <- Dead_Table %>% filter(BC7 != "2632731" & BC7 != "2632726" & BC7 != "2632716" & BC7 != "2632721")
Dead_Table <- rbind(Dead_Table, CT_Day7)

Dead_Table["0"][is.na(Dead_Table["0"])] <- 0
Dead_Table$`7` <- ifelse(is.na(Dead_Table$`7`), Dead_Table$Total_Flies, Dead_Table$`7`)
Dead_Table$`6` <- ifelse(is.na(Dead_Table$`6`), round((Dead_Table$`7`+Dead_Table$`5`)/2), Dead_Table$`6`)
Dead_Table$`6` <- ifelse(is.na(Dead_Table$`6`), Dead_Table$Total_Flies, Dead_Table$`6`)
Dead_Table$`5` <- ifelse(is.na(Dead_Table$`5`), round((Dead_Table$`6`+Dead_Table$`4`)/2), Dead_Table$`5`)
Dead_Table$`5` <- ifelse(is.na(Dead_Table$`5`), Dead_Table$Total_Flies, Dead_Table$`5`)
Dead_Table$`4` <- ifelse(is.na(Dead_Table$`4`), round((Dead_Table$`5`+Dead_Table$`3`)/2), Dead_Table$`4`)
Dead_Table$`4` <- ifelse(is.na(Dead_Table$`4`), Dead_Table$Total_Flies, Dead_Table$`4`)

#Create Change in Dead
Change_Dead_Table <- Dead_Table[c(1:5)]
Change_Dead_Table$'1' <- Dead_Table$`1`
Change_Dead_Table$'2' <- Dead_Table$`2`-Dead_Table$`1`
Change_Dead_Table$'3' <- Dead_Table$`3`-Dead_Table$`2`
Change_Dead_Table$'4' <- Dead_Table$`4`-Dead_Table$`3`
Change_Dead_Table$'5' <- Dead_Table$`5`-Dead_Table$`4`
Change_Dead_Table$'6' <- Dead_Table$`6`-Dead_Table$`5`
Change_Dead_Table$'7' <- Dead_Table$`7`-Dead_Table$`6`

Change_Dead_Long <- Change_Dead_Table %>% pivot_longer(c("1":"7"),names_to = "Day", values_to = "Change_Dead")
head(Change_Dead_Long)
#   BC7     Replicate Population Treatment SetUpBlockNumeric Day   Change_Dead
#   <chr>   <chr>     <chr>      <chr>     <chr>             <chr>       <dbl>
# 1 2632301 B10       SelZ       TrtZ      1                 1               0
# 2 2632301 B10       SelZ       TrtZ      1                 2               4
# 3 2632301 B10       SelZ       TrtZ      1                 3              13
# 4 2632301 B10       SelZ       TrtZ      1                 4               2
# 5 2632301 B10       SelZ       TrtZ      1                 5               1
# 6 2632301 B10       SelZ       TrtZ      1                 6               0

#Combine counts for each day based on Replicate, Population, Treatment, and generation (SetupBlockNumeric)
Change_Dead_Pop <- aggregate(Change_Dead~Replicate+Population+Treatment+Day+SetUpBlockNumeric, Change_Dead_Long, sum)
head(Change_Dead_Pop)
#   Replicate Population Treatment Day SetUpBlockNumeric Change_Dead
# 1       B10       SelC      TrtC   1                 1           0
# 2        B9       SelC      TrtC   1                 1           0
# 3       B10       SelZ      TrtC   1                 1           0
# 4        B9       SelZ      TrtC   1                 1           0
# 5       B10       SelC      TrtZ   1                 1           3
# 6        B9       SelC      TrtZ   1                 1           4


#Time of Death for each fly####
#Create a table where each fly's death becomes a data point
#Create a TOD, Population, Treatment, Replicate, and Generation vector
#Code is doing: Take the value in the Day's column, and repeat x times, with
#x being the value in the Change_Dead column.
TOD_Rep <- rep(Change_Dead_Pop$Day, Change_Dead_Pop$Change_Dead)
#Doing same as below except switching out Day's column. 
Pop_Rep <- rep(Change_Dead_Pop$Population, Change_Dead_Pop$Change_Dead)
Treat_Rep <- rep(Change_Dead_Pop$Treatment, Change_Dead_Pop$Change_Dead)
Rep_Rep <- rep(Change_Dead_Pop$Replicate, Change_Dead_Pop$Change_Dead)
Gen_Rep <- rep(Change_Dead_Pop$SetUpBlockNumeric, Change_Dead_Pop$Change_Dead)
#Because we are repeating by the number in the Change_Dead column, when we compare
#the head() of Change_Dead_Pop and TOD, Change_Dead_Pop first 4 rows have the 
#control treat, but TOD begins with zinc treatment. 

#cbind vectors together
TOD <- cbind(Pop_Rep, Treat_Rep)
TOD <- cbind(TOD, Rep_Rep)
TOD <- cbind(TOD,Gen_Rep)
TOD <- cbind(TOD, TOD_Rep)
TOD <- as.data.frame(TOD)
TOD$TOD_Rep <- as.numeric(TOD$TOD_Rep)
colnames(TOD) <- c("Selection","Treatment","Replicate","Gen","TOD")
head(TOD)
#   Selection Treatment Replicate Gen TOD
# 1      SelC      TrtZ       B10   1   1
# 2      SelC      TrtZ       B10   1   1
# 3      SelC      TrtZ       B10   1   1
# 4      SelC      TrtZ        B9   1   1
# 5      SelC      TrtZ        B9   1   1
# 6      SelC      TrtZ        B9   1   1

#Give more informative Gen values
TOD$Gen[TOD$Gen == "1"] <- c("Gen. 1")
TOD$Gen[TOD$Gen == "2"] <- c("Gen. 3")

#Check for control treatment death####
TOD %>% filter(Treatment == "TrtC")

#  Selection Treatment Replicate    Gen TOD
#1      SelC      TrtC        B9 Gen. 1   7
#Out of 310 flies tested on control media only one died at Day 7. 
#Each population had 155 flies tested on control media. 
Control_Treatment <- Dead_Table %>% filter(Treatment == "TrtC")
Control_Treatment_Combine <- aggregate(Total_Flies~Population, Control_Treatment, sum)
head(Control_Treatment_Combine)
#   Population Total_Flies
# 1       SelC        155
# 2       SelZ        155
SelC_Control <- rep(0, 155)
SelZ_Control <- rep(0, 155)
SelC_Control[1] <- 1
t.test(SelC_Control, SelZ_Control)
# Welch Two Sample t-test
# 
# data:  SelC_Control and SelZ_Control
# t = 1, df = 154, p-value = 0.3189
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.006293471  0.019196697
# sample estimates:
#   mean of x   mean of y 
# 0.006451613 0.000000000 
#There is not significant difference between the zinc selected and non-selected 
#flies on control media. 

#Due to only 1 death, this death happening at the end of the experiment, and no 
#significant differences between the populations, will do TOD analysis without
#this death. 

#TOD Analysis####
TOD_TrtZ <- TOD %>% filter(Treatment == "TrtZ")
colnames(TOD_TrtZ)[1] <- c("Population")

TOD_aov <- aov(lm(TOD~Gen*Replicate*Population, data = TOD_TrtZ))
anova(TOD_aov)
# Response: TOD
#                            Df  Sum Sq Mean Sq F value    Pr(>F)    
# Gen                         1   37.36  37.357 39.2059 4.509e-10 ***
# Replicate                   1   27.50  27.499 28.8593 8.538e-08 ***
# Population                  1    6.43   6.431  6.7491  0.009437 ** 
# Gen:Replicate               1    9.89   9.893 10.3830  0.001289 ** 
# Gen:Population              1    1.81   1.812  1.9019  0.167994    
# Replicate:Population        1    0.01   0.009  0.0089  0.924667    
# Gen:Replicate:Population    1    0.21   0.212  0.2222  0.637401    
# Residuals                2384 2271.60   0.953      

#Generation 1 alone analysis
Gen1_TOD <- TOD_TrtZ %>% filter(Gen == "Gen. 1")
G1_TOD_aov <- aov(lm(TOD~Replicate*Population, data = Gen1_TOD))
anova(G1_TOD_aov)
# Response: TOD
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Replicate              1   0.21 0.20774  0.1948 0.6591
# Population             1   0.13 0.13196  0.1237 0.7251
# Replicate:Population   1   0.10 0.10408  0.0976 0.7548
# Residuals            792 844.62 1.06643

#Generation 3 alone analysis
Gen3_TOD <- TOD_TrtZ %>% filter(Gen == "Gen. 3")
G3_TOD_aov <- aov(lm(TOD~Replicate*Population, data = Gen3_TOD))
anova(G3_TOD_aov)
# Response: TOD
#                        Df  Sum Sq Mean Sq F value   Pr(>F)    
# Replicate               1   37.18  37.184 41.4842 1.57e-10 ***
# Population              1    8.11   8.111  9.0493 0.002669 ** 
# Replicate:Population    1    0.12   0.116  0.1296 0.718884    
# Residuals            1592 1426.99   0.896 

#Save TOD data table
write.csv(TOD, "X-QTL Subpopulations TOD.csv")
