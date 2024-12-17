#Introduction####
#Title Supplementary Figure 5
#Create: 12/13/24
#Last Updated: 12/17/24, KMH
#Packages:
library(tidyverse)
library(Rmisc)

#For this analysis will need "Phenotyping Cohorts TOD.csv" and "X-QTL Subpopulations TOD.csv"

#Read in Data####
cohorts_TOD <- read.csv("Phenotyping Cohorts TOD.csv", header = TRUE)
subpopulations_TOD <- read.csv("X-QTL Subpopulations TOD.csv", header = TRUE)

#Adjust names and columns so have consistent nomenclature and same number columns####
subpopulations_TOD <- subpopulations_TOD[-c(1)]
subpopulations_TOD$Gen[subpopulations_TOD$Gen == "Gen. 1"] <- "X-QTL Subpopulations Gen. 1" 
subpopulations_TOD$Gen[subpopulations_TOD$Gen == "Gen. 3"] <- "X-QTL Subpopulations Gen. 3" 
subpopulations_TOD$Selection[subpopulations_TOD$Selection == "SelZ"] <- "Zinc Selected"
subpopulations_TOD$Selection[subpopulations_TOD$Selection == "SelC"] <- "Control Non-Selected"
subpopulations_TOD <- subpopulations_TOD %>% filter(Treatment == "TrtZ")
subpopulations_TOD$Replicate <- gsub("B","R", subpopulations_TOD$Replicate)
colnames(subpopulations_TOD)[1] <- "Population"

cohorts_TOD <- cohorts_TOD[-c(1)]
cohorts_TOD$Population[cohorts_TOD$Population == "Control"] <- "Control Non-Selected"
cohorts_TOD$Population[cohorts_TOD$Population == "Zinc"] <- "Zinc Selected"
cohorts_TOD$Replicate <- "R1"
cohorts_TOD$Gen <- "Phenotyping-Only Cohort"

TOD_all <- rbind(subpopulations_TOD, cohorts_TOD)

#All populations were tested on 100mM ZnCl2

#Statiscal Test Results####
phenotyping_cohort_ttest <- t.test(cohorts_TOD %>% filter(Population == "Zinc Selected") %>% pull(TOD),
                                   cohorts_TOD %>% filter(Population != "Zinc Selected") %>% pull(TOD))

subpop_Gen1_R9 <- t.test(subpopulations_TOD %>% filter(Population =="Zinc Selected" & Replicate == "R9" & Gen == "X-QTL Subpopulations Gen. 1") %>% pull(TOD),
                         subpopulations_TOD %>% filter(Population !="Zinc Selected" & Replicate == "R9" & Gen == "X-QTL Subpopulations Gen. 1") %>% pull(TOD))

subpop_Gen1_R10 <- t.test(subpopulations_TOD %>% filter(Population =="Zinc Selected" & Replicate == "R10" & Gen == "X-QTL Subpopulations Gen. 1") %>% pull(TOD),
                          subpopulations_TOD %>% filter(Population !="Zinc Selected" & Replicate == "R10" & Gen == "X-QTL Subpopulations Gen. 1") %>% pull(TOD))

subpop_Gen3_R9 <- t.test(subpopulations_TOD %>% filter(Population =="Zinc Selected" & Replicate == "R9" & Gen == "X-QTL Subpopulations Gen. 3") %>% pull(TOD),
                         subpopulations_TOD %>% filter(Population !="Zinc Selected" & Replicate == "R9" & Gen == "X-QTL Subpopulations Gen. 3") %>% pull(TOD))

subpop_Gen3_R10 <- t.test(subpopulations_TOD %>% filter(Population =="Zinc Selected" & Replicate == "R10" & Gen == "X-QTL Subpopulations Gen. 3") %>% pull(TOD),
                          subpopulations_TOD %>% filter(Population !="Zinc Selected" & Replicate == "R10" & Gen == "X-QTL Subpopulations Gen. 3") %>% pull(TOD))

#Plot####
label_stats <- as.data.frame(rbind(cbind("Phenotyping-Only Cohort", "7.708e-11", 4.1, "R1",""),
                                  cbind("X-QTL Subpopulations Gen. 1", round(subpop_Gen1_R9[["p.value"]],3), 4.1, "R9",""),
                                  cbind("X-QTL Subpopulations Gen. 1", round(subpop_Gen1_R10[["p.value"]],3), 4.1, "R10",""),
                                  cbind("X-QTL Subpopulations Gen. 3", round(subpop_Gen3_R9[["p.value"]],3), 4.1, "R9",""),
                                  cbind("X-QTL Subpopulations Gen. 3", round(subpop_Gen3_R10[["p.value"]],3), 4.1, "R10","")))
colnames(label_stats) <- c("Gen", "P-value", "TOD", "Replicate","Population")
label_stats$TOD <- as.numeric(label_stats$TOD)
label_stats$Gen <- factor(label_stats$Gen, levels = c("X-QTL Subpopulations Gen. 1",
                                                    "X-QTL Subpopulations Gen. 3",
                                                    "Phenotyping-Only Cohort"))

TOD_all_ci <- summarySE(data=TOD_all, measurevar = "TOD", groupvars = c("Population","Replicate","Gen"))

TOD_all_ci$Population <- factor(TOD_all_ci$Population,levels = c("Control Non-Selected","Zinc Selected",""))
label_stats$Population <- factor(label_stats$Population, levels = c("Control Non-Selected","Zinc Selected",""))
TOD_all_ci$Gen <- factor(TOD_all_ci$Gen, levels = c("X-QTL Subpopulations Gen. 1",
                                                     "X-QTL Subpopulations Gen. 3",
                                                     "Phenotyping-Only Cohort"))
TOD_all_ci$Replicate <- factor(TOD_all_ci$Replicate, levels = c("R1","R9","R10"))

TOD_adults <- 
  TOD_all_ci %>% ggplot(aes(x=Replicate, y=TOD, color=Population))+
  geom_errorbar(aes(ymin=TOD-ci, ymax=TOD+ci), width=0.3,linewidth=1, position = position_dodge(.3))+
  geom_point(size=3, position = position_dodge(.3))+theme_bw()+
  scale_color_manual(values = c("#00AFBB","#FC4E07","Black"))+
  #ggtitle(title)+
  facet_wrap(~Gen, scales = "free_x")+
  geom_text(data=label_stats,label=paste("T-test,",label_stats$`P-value`, sep = " "), show.legend = FALSE, color="Black")+
  ylab("Time of Death (Days)")+
  theme(legend.position = "bottom", text = element_text(size = 10))

ggsave("Phenotyping_Population_Adults.png",TOD_adults,width = 7, height = 4, units = "in")

