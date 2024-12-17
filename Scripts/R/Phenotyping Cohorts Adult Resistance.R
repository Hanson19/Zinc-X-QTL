#Introduction####
#Title: X-QTL Phenotyping Only Cohort Adult Resistance
#Purpose: To analyze time of death of adult zinc naive female flies derived
#from zinc selected and non-selected flies  generated in same design as x-QTL design
#Created: 11/28/22
#Last Edited: 7/05/23
#Packages:
library(tidyverse)
library(Rmisc)

#For this analysis you will need PP-II_G2_Adult_Barcodes_11182022 and 
#"PPII_G2_Allscans.txt"

#Read in Data####
AllScans <- read.table("PPII_G2_Allscans.txt", header = TRUE)
head(AllScans)
#       BC7 NumDead DateDDMMYY TimeHHMM ScanYYMMDDHHMM ScanExperimenter
# 1 2800350       0     191122     1106     2211191106              022
# 2 2800369       0     191122     1106     2211191106              022
# 3 2800304       0     191122     1106     2211191106              022
# 4 2800354       0     191122     1106     2211191106              022
# 5 2800398       0     191122     1106     2211191106              022
# 6 2800348       0     191122     1106     2211191106              022

#AllScans has every vial scan that counted the number of dead.
#Each unique BC7 represents 1 vial.
#NumDead is the number of dead flies in vial at time of scanning
#DateDDMMYY is the day the vials were setup
#TimeHHMM is the time that vials were first scanned during setup. 
#ScanYYMMDDHHMM is the date and time of the scan to count the dead was done.
#22111941106 = November 11th, 2022 and 11:06am.
#ScanExperimenter is a remnant column from previous code (See PP-1 Adult Analysis)
#All scans were done by the same experimenter and the 022 comes from the last 
#three characters of the name of the scans which stands for the date. 

#Check for Inconsistencies in Scans####
Barcodes <- unique(AllScans$BC7)
AllScans <- AllScans %>% arrange(DateDDMMYY)

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

AllScans <- findinconsistency.func(AllScans,Barcodes)

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
  ScansBy7BC <- AllScans %>% filter(BC7 == i)
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

dim(AllScans)
#504 7
dim(Gold)
#504 8
dim(Valid)
#NULL
dim(Manual)
#NULL
#There are no inconsistencies between scans.
#Because there are none will continue using AllScans

#Correct for Initial Dead####
#Convert date numbers to dates and then factor them
AllScans$Formatted_Date <- strptime(AllScans$DateDDMMYY, "%d%m%y")
AllScans$Formatted_Date <- factor(AllScans$Formatted_Date, 
                                  levels = c("2022-11-18", "2022-11-19", 
                                             "2022-11-20", "2022-11-21",
                                             "2022-11-22", "2022-11-23", 
                                             "2022-11-24", "2022-11-25",
                                             "2022-11-26"))
AllScans$NumDead <- as.numeric(AllScans$NumDead)

#Format into wide table so each vial is its own row
Total_Dead_Table <- AllScans %>% pivot_wider(id_cols = BC7, 
                                             names_from = Formatted_Date, 
                                             values_from = NumDead)

#Rename 2022-11-18 to Initial Dead and 2022-11-26 to Total Dead
colnames(Total_Dead_Table)[2] <- c("Initial_Dead")
colnames(Total_Dead_Table)[10] <- c("Total_Flies")
head(Total_Dead_Table)
#       BC7 Initial_Dead `2022-11-19` `2022-11-20` `2022-11-21` `2022-11-22` `2022-11-23` `2022-11-24` `2022-11-25` Total_Flies
#     <int>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>       <dbl>
# 1 2800349            0            0            0           10           18           20           20           20          20
# 2 2800305            0            0            1            7           16           20           20           20          20
# 3 2800353            0            0            4           11           19           20           20           20          20
# 4 2800399            0            0            2           16           19           20           20           20          20
# 5 2800308            0            0            0            7           17           19           20           20          20
# 6 2800303            0            0            1            8           18           20           20           20          20

#Correct for Initial Dead
Total_Dead_Table %>% filter(BC7 == "2800310")
#         BC7 Initial_Dead `2022-11-19` `2022-11-20` `2022-11-21` `2022-11-22` `2022-11-23` `2022-11-24` `2022-11-25` Total_Flies
#       <int>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>       <dbl>
#   1 2800310            1            1            3            8           18           20           20           20          20

#Subtract the value in Initial_Dead column from all other columns
Total_Dead_Table$`2022-11-19` <- Total_Dead_Table$`2022-11-19`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-20` <- Total_Dead_Table$`2022-11-20`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-21` <- Total_Dead_Table$`2022-11-21`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-22` <- Total_Dead_Table$`2022-11-22`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-23` <- Total_Dead_Table$`2022-11-23`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-24` <- Total_Dead_Table$`2022-11-24`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$`2022-11-25` <- Total_Dead_Table$`2022-11-25`-Total_Dead_Table$Initial_Dead
Total_Dead_Table$Total_Flies <- Total_Dead_Table$Total_Flies-Total_Dead_Table$Initial_Dead

Total_Dead_Table %>% filter(BC7 == "2800310")
#         BC7 Initial_Dead `2022-11-19` `2022-11-20` `2022-11-21` `2022-11-22` `2022-11-23` `2022-11-24` `2022-11-25` Total_Flies
#       <int>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>        <dbl>       <dbl>
#   1 2800310            1            0            2            7           17           19           19           19          19

#Get Barcode ID information in table
#Read in Barcode Information
Total_Dead_Table$BC7 <- as.character(Total_Dead_Table$BC7)
Barcodes <- read.table("PP-II G2 Adult Barcodes 11182022.txt", header = TRUE)
colnames(Barcodes)[1] <- c("BC7")
Barcodes$BC7 <- as.character(Barcodes$BC7)
Total_Dead_Table <- right_join(Barcodes, Total_Dead_Table,by=c("BC7"))
head(Total_Dead_Table)
#       BC7 Population  Treatment Setup_Date Initial_Dead 2022-11-19 2022-11-20 2022-11-21 2022-11-22 2022-11-23 2022-11-24 2022-11-25 Total_Flies
# 1 2800350    Control        H2O 2022-11-18            0          0          0          0          0          0          0          0          19
# 2 2800304       Zinc        H2O 2022-11-18            0          0          0          0          0          0          0          0          19
# 3 2800345    Control 100mM_Zinc 2022-11-18            0          0          0         12         18         20         20         20          20
# 4 2800309       Zinc 100mM_Zinc 2022-11-18            0          0          1          6         18         20         20         20          20
# 5 2800354    Control        H2O 2022-11-18            0          0          0          0          0          0          0          0          20
# 6 2800398       Zinc        H2O 2022-11-18            0          0          0          0          0          0          0          0          20

#Now have identifying information about each vial.
#Population: Are the flies derived from the control non-selected population 
#(Control) or zinc selected population (zinc)
#Treatment: Were flies put on control H2O media or media supplemented with
#100mM ZnCl2
#Setup_Date: Date flies were put into vials. All flies were put in on the same
#day 2022-11-18

#Change in Dead####
#Pull out identifying information and first day of counting
Change_Dead <- Total_Dead_Table[,c(1:4,6)]
#Subtract previous day's count 
Change_Dead$'2022-11-20' <- Total_Dead_Table$`2022-11-20`-Total_Dead_Table$`2022-11-19`
Change_Dead$'2022-11-21' <- Total_Dead_Table$`2022-11-21`-Total_Dead_Table$`2022-11-20`
Change_Dead$'2022-11-22' <- Total_Dead_Table$`2022-11-22`-Total_Dead_Table$`2022-11-21`
Change_Dead$'2022-11-23' <- Total_Dead_Table$`2022-11-23`-Total_Dead_Table$`2022-11-22`
Change_Dead$'2022-11-24' <- Total_Dead_Table$`2022-11-24`-Total_Dead_Table$`2022-11-23`
Change_Dead$'2022-11-25' <- Total_Dead_Table$`2022-11-25`-Total_Dead_Table$`2022-11-24`

#Create Long Form Change of Dead table and then sum based on population and treatment and Date 
Change_Dead_Long <- Change_Dead %>% pivot_longer(c(`2022-11-19`:`2022-11-25`), 
                                                 names_to = "Date", values_to = "Chang_Dead")
dim(Change_Dead_Long)
#[1] 392 6

#Sum all the counts based on population,treatment and date
Change_Dead_Long_Pop <- aggregate(Chang_Dead~Population+Treatment+Setup_Date+Date, Change_Dead_Long, sum)                                                 
dim(Change_Dead_Long_Pop)
#[1] 28 5

#Time of Death for each fly####
#Get the difference in time (Days) between setup date and scan date
Change_Dead_Long_Pop$TOD <- as.numeric(difftime(Change_Dead_Long_Pop$Date, Change_Dead_Long_Pop$Setup_Date))
#   Population  Treatment Setup_Date       Date Chang_Dead TOD
# 1    Control 100mM_Zinc 2022-11-18 2022-11-19          0   1
# 2       Zinc 100mM_Zinc 2022-11-18 2022-11-19          0   1
# 3    Control        H2O 2022-11-18 2022-11-19          0   1
# 4       Zinc        H2O 2022-11-18 2022-11-19          0   1
# 5    Control 100mM_Zinc 2022-11-18 2022-11-20         21   2
# 6       Zinc 100mM_Zinc 2022-11-18 2022-11-20         41   2

#Create a table where each fly's death becomes a data point
#Create a TOD, Population, Treatment, Replicate, and Generation vector
#Code is doing: Take the value in the TOD column, and repeat x times, with
#x being the value in the Chang_Dead column.
TOD_Rep <- rep(Change_Dead_Long_Pop$TOD, Change_Dead_Long_Pop$Chang_Dead)
#Doing same as below except switching out Day's column. 
Pop_Rep <- rep(Change_Dead_Long_Pop$Population, Change_Dead_Long_Pop$Chang_Dead)
Treat_Rep <- rep(Change_Dead_Long_Pop$Treatment, Change_Dead_Long_Pop$Chang_Dead)

#cbind vectors together
TOD <- cbind(Pop_Rep, Treat_Rep)
TOD <- cbind(TOD, TOD_Rep)
TOD <- as.data.frame(TOD)
TOD$TOD_Rep <- as.numeric(TOD$TOD_Rep)
colnames(TOD) <- c("Population","Treatment","TOD")
head(TOD)
#   Population  Treatment TOD
# 1    Control 100mM_Zinc   2
# 2    Control 100mM_Zinc   2
# 3    Control 100mM_Zinc   2
# 4    Control 100mM_Zinc   2
# 5    Control 100mM_Zinc   2
# 6    Control 100mM_Zinc   2
#Each row represent a fly and has identifying information about this fly such
#as what population it came from (control non-selected or zinc selected), what
#treatment were the flies on (Only flies on 100mM ZnCl2 died during this experiment)
#so all are "100mM_Zinc", and the the flies time of death (TOD) which how many 
#days post setup did the fly die. 


#TOD Analyze
#Only variable to consider is population so can just do a t-test.
Control <- TOD %>% filter(Population == "Control")
Zinc <- TOD %>% filter(Population == "Zinc")
t.test(Control$TOD, Zinc$TOD)
# Welch Two Sample t-test
# 
# data:  Control$TOD and Zinc$TOD
# t = 6.5778, df = 993.75, p-value = 7.708e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 0.2607484 0.4824748
# sample estimates:
# mean of x mean of y 
# 3.936128  3.564516 

#Zinc selected population dies faster on 100mM ZnCl2 then the control non-selected
#population 

write.csv(TOD, "Phenotyping Cohorts TOD.csv")
