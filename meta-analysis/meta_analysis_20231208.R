### META-ANALYSIS OF PATHWAY-SPECIFIC PRS ###

#---Load libaries---####
library(tidyverse)
library(data.table)
library(meta)
library(grid)


#---Read in data---####

nih_cs <- read_csv("./NIH/report/cs.csv")

nih_lt <- read_csv("./NIH/report/lt.csv")

nih_surv <- read_csv("./NIH/report/surv.csv")

amp_agediagnosis <- read.table("./AMP/pathways_agediagnosis_results.txt", header = T)

amp_dementia <- read.table("./AMP/pathways_dementia_results.txt", header = T)

amp_hy3 <- read.table("./AMP/pathways_hy3_results.txt", header = T)

amp_mdsupdrs2 <- read.table("./AMP/pathways_mdsupdrs2_results.txt", header = T)

amp_mdsupdrs3 <- read.table("./AMP/pathways_mdsupdrs3_results.txt", header = T)

amp_moca <- read.table("./AMP/pathways_moca_results.txt", header = T)

amp_PDQ8 <- read.table("./AMP/pathways_PDQ8_results.txt", header = T)

amp_RBD <- read.table("./AMP/pathways_RBD_results.txt", header = T)

Oslo_aao <- read.table("./Oslo/aao_pathways_results_PRSChang0.05.txt", header = T)

Oslo_hy3 <- read.table("./Oslo/HY3_pathways_results_PRSChang0.05.txt", header = T)

Oslo_mortality <- read.table("./Oslo/mortality_pathways_results_PRSChang0.05.txt", header = T)

Oslo_mdsupdrs3 <- read.table("./Oslo/updrs3_pathways_results_PRSChang0.05.txt", header = T)

McGill_aao <- read.table("./McGill/aao_pathways_results_PRSChang0.05.txt", header = T)

McGill_hy3 <- read.table("./McGill/HY3_pathways_results_PRSChang0.05.txt", header = T)

McGill_RBD <- read.table("./McGill/RBD_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_aao <- read.table("./PROBAND/aao_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_dementia <- read.table("./PROBAND/dementia_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_hy3 <- read.table("./PROBAND/HY3_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_mdsupdrs3 <- read.table("./PROBAND/UPDRSIII_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_mdsupdrs2 <- read.table("./PROBAND/UPDRSII_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_moca <- read.table("./PROBAND/MOCA_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_RBD <- read.table("./PROBAND/RBD_pathways_results_PRSChang0.05.txt", header = T)

PROBAND_PDQ8 <- read.table("./PROBAND/PDQ8_pathways_results_PRSChang0.05.txt", header = T)

#---Read in Fox Insight results---####

fox_aao <- read.table("./FoxInsight/pathways_results_PRSChang0.05_aao.txt", header = T)

fox_mdsupdrs2 <- read.table("./FoxInsight/pathways_results_PRSChang0.05_updrs2.txt", header = T)

#---Make function to run meta-analysis and plot forest plots---####

#Make function to run random effects meta-analysis and forest plot on specific outcome
#Loops over all pathways
#The outcome must be named as the variables are named in the datasets
run_meta <- function(outcome, data) {
  
  #Make results dataframe for all pathways
  meta_results <- as.data.frame(matrix(ncol = 9))
  colnames(meta_results) <- c("pathway", "beta", "se", "95CI", "pval", "n_studies", "I2", "CochransQ_pval", "N_inds")
  
  #Make list of pathways
  pathways <- c("adaptive_immune", "alpha_synuclein", "innate_immune", "lysosome", "endocytosis",
                "microglia", "monocytes", "mitochondria")
  
  for (nm in 1:length(pathways)) {
    
    #Filter the merged dataset for this pathway
    #Remove cohorts missing data
    filtered <- data %>% 
      filter(pathway ==  pathways[nm]) %>% 
      drop_na()
  
    #Run meta-analysis using random-effects model
    meta <- metagen(Coeff,
                    se,
                    data = filtered,
                    studlab = dataset,
                    fixed = TRUE,
                    random = TRUE,
                    prediction = TRUE,
                    sm = "SMD") #all effect sizes are in betas
    
    #Put results into results dataframe for each pathway
    meta_results[nm,1] <- pathways[nm] #name of pathway
    meta_results[nm,2] <- meta$TE.random #random effect effect size
    meta_results[nm,3] <- meta$seTE.random #random effect effect size SE
    meta_results[nm,4] <- paste(sprintf("%.2f", round(meta$lower.random, 2)), ", ", sprintf("%.2f", round(meta$upper.random,2)), sep = "") #lower to upper 95% CI
    meta_results[nm,5] <- meta$pval.random #pvalue 
    meta_results[nm,6] <- meta$k.study #number of studies
    meta_results[nm,7] <- meta$I2 #Isquared
    meta_results[nm,8] <- meta$pval.Q #Cochran's Q p-value
    meta_results[nm,9] <- sum(filtered$N)  #Calculate total number of individuals
    
    #Make forest plot
    filename_plot <- paste("./plots/forestplot_", outcome, "_", pathways[nm], ".png", sep = "")
    png(file = filename_plot, width = 1000, height = 700) 
    forest.meta(meta,
           layout = "meta",
           leftcols = c("studlab"),
           rightcols = c("effect.ci"),
           rightlabs = c("SMD", "[95% CI]"),
           colgap.forest.left = unit(25,"mm"),
           print.tau2 = FALSE,
           prediction = FALSE,
           col.diamond = "blue",
           fontsize = 18,
           spacing = 1.7,
           squaresize = 0.9,
           lwd = 2,
           sortvar = studlab
    )
    text <- paste("Effect size for ", pathways[nm], " PRS on ", outcome, sep = "")
    grid.text(text, 0.5, .98, 
              gp=gpar(fontsize = 8, cex=2, fontface = 2))
    dev.off()
    
    
  }
  
  #Remove empty rows
  data_nomiss <- data %>% 
    na.omit()
  
  #Save raw cohort results
  write.table(data_nomiss, paste("./outputs/cohorts_", outcome, ".txt", sep = ""),
                quote = F, row.names = F, col.names = T, sep = "\t")
  
  return(meta_results)
}

#---Meta-analysis MDS-UPDRS3---####

nih_mdsupdrs3 <- nih_lt %>% 
  filter(outcome == "MDS_UPDRS3") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

amp_mdsupdrs3 <- amp_mdsupdrs3 %>%
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N = N_inds)

Oslo_mdsupdrs3 <- Oslo_mdsupdrs3 %>% 
  mutate(dataset = "Oslo") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

PROBAND_mdsupdrs3 <- PROBAND_mdsupdrs3 %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)


mdsupdrs3 <- rbind(nih_mdsupdrs3, amp_mdsupdrs3, Oslo_mdsupdrs3, PROBAND_mdsupdrs3)

meta_mdsupdrs3 <- run_meta("mdsupdrs3", mdsupdrs3)


#---Meta-analysis age at diagnosis---####

nih_aao <- nih_cs %>% 
  filter(outcome == "AAO") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

amp_aao <- amp_agediagnosis %>% 
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N)

Oslo_aao <- Oslo_aao %>% 
  mutate(dataset = "Oslo")

McGill_aao <- McGill_aao %>% 
  mutate(dataset = "QPN")

PROBAND_aao <- PROBAND_aao %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

Fox_aao <- Fox_aao %>% 
  mutate(dataset = "Fox") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

aao <- rbind(nih_aao, amp_aao, Oslo_aao, McGill_aao, PROBAND_aao)

meta_aao <- run_meta("AAO", aao)


#---Meta-analysis HY3+---####

nih_hy3 <- nih_surv %>% 
  filter(outcome == "HY3") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

amp_hy3 <- amp_hy3 %>% 
  select(dataset = study, pathway, Coeff, se, Pvalue, N)

Oslo_hy3 <- Oslo_hy3 %>% 
  mutate(dataset = "Oslo") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

McGill_hy3 <- McGill_hy3 %>% 
  mutate(dataset = "QPN") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)
  
PROBAND_hy3 <- PROBAND_hy3 %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

hy3 <- rbind(nih_hy3, amp_hy3, Oslo_hy3, McGill_hy3, PROBAND_hy3)

#Remove duplicates - data for PRECEPT is duplicated for endocytosis pathway
hy3_unique <- hy3 %>% 
  distinct(dataset, pathway, .keep_all = TRUE)

meta_hy3 <- run_meta("hy3", hy3_unique)

#---Meta-analysis dementia---####

nih_dementia <- nih_surv %>% 
  filter(outcome == "MCI") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

amp_dementia <- amp_dementia %>% 
  select(dataset = study, pathway, Coeff, se, Pvalue, N)

PROBAND_dementia <- PROBAND_dementia %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

dementia <- rbind(nih_dementia, amp_dementia, PROBAND_dementia)

meta_dementia <- run_meta("dementia", dementia)

#---Meta-analysis MDS-UPDRS2---####

nih_mdsupdrs2 <- nih_lt %>% 
  filter(outcome == "MDS_UPDRS2") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

amp_mdsupdrs2 <- amp_mdsupdrs2 %>% 
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N = N_inds)

PROBAND_mdsupdrs2 <- PROBAND_mdsupdrs2 %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

fox_mdsupdrs2 <- fox_mdsupdrs2 %>% 
  mutate(dataset = "Fox") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)
  
mdsupdrs2 <- rbind(nih_mdsupdrs2, amp_mdsupdrs2, PROBAND_mdsupdrs2)

#Run analysis
meta_mdsupdrs2 <- run_meta("mdsupdrs2", mdsupdrs2)

#---Meta-analysis MOCA---####

nih_moca <- nih_lt %>% 
  filter(outcome == "MOCA") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

amp_moca <- amp_moca %>% 
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N = N_inds)

PROBAND_moca <- PROBAND_moca %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

moca <- rbind(nih_moca, amp_moca, PROBAND_moca)

#Run analysis 
meta_moca <- run_meta("MOCA", moca)

#---Meta-analysis PDQ8---####

amp_PDQ8 <- amp_PDQ8 %>% 
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N = N_inds)

PROBAND_PDQ8 <- PROBAND_PDQ8 %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

PDQ8 <- rbind(amp_PDQ8, PROBAND_PDQ8)

#Run analysis 
meta_PDQ8 <- run_meta("PDQ8", PDQ8)


#---Meta-analysis RBD---####

nih_rbd <- nih_cs %>% 
  filter(outcome == "pRBD") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

amp_rbd <- amp_RBD %>% 
  select(dataset = cohort, pathway, Coeff, se, Pvalue, N = N_inds)

McGill_rbd <- McGill_RBD %>% 
  mutate(dataset = "QPN") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N)

PROBAND_rbd <- PROBAND_RBD %>% 
  mutate(dataset = "PROBAND") %>% 
  select(dataset, pathway, Coeff, se, Pvalue, N = N_inds)

rbd <- rbind(nih_rbd, amp_rbd, McGill_rbd, PROBAND_rbd)

#Run analysis 
meta_rbd <- run_meta("RBD", rbd)


#---Write results files---####

list_results <- grep("meta_",names(.GlobalEnv),value=TRUE)
Pattern1_list <- do.call("list",mget(list_results))


for (i in 1:length(list_results)) {
  
  #Get name of dataframe
  results <- list_results[i]

  today <- Sys.Date()
  
  #Filename
  filename <- paste("./outputs/", results, "_", today, ".txt", sep = "")
  
  #Write as output
  write.table(get(results), filename, quote = F, row.names = F, col.names = T, sep = "\t")
}



##### SIGNIFICANCE FDR #####

#Combine all pathways 

outcomes <- c("aao", "dementia", "hy3", 'mdsupdrs2', "mdsupdrs3", "moca", "PDQ8", "rbd")

datalist <- list()

for (i in 1:length(outcomes)) {
  dataframe_name <- paste("meta_", outcomes[i], sep = "")
  df <- get(dataframe_name)
  
  df_selected <- df %>%
    mutate(outcome = outcomes[i]) %>% 
    select(outcome, pathway, RE_SMD, pval)
  
  datalist[[i]] <- df_selected
}


#Combine
combined_results <- do.call(rbind, datalist)

#Sort by p-value
combined_results_sorted <- combined_results %>% 
  arrange(pval)

#Make rank column
combined_results_sorted$rank <- rank(combined_results_sorted$pval)

#Make significance column - using Benjamini & Hochberg FDR
combined_results_sorted <- combined_results_sorted %>% 
  mutate(Pvalue_adjusted = rank/64 * 0.05) %>% 
  mutate(significant = ifelse(pval < Pvalue_adjusted, "significant", "NS"))


##### META-ANALYSE WITH FOX INSIGHT #####

#---Beta-beta-plots---####

Fox_aao <- read.table("../../FoxDEN/outputs/pathways_results_PRSChang0.05_aao.txt", header = T)

#Add column for dataset
Fox_aao <- Fox_aao %>% 
  mutate(dataset = "FoxInsight")

#Combine with clinical cohorts results
aao_Fox <- rbind(aao, foxinsight_results)

#Filter the merged dataset for this pathway
#Remove cohorts missing data
aao_Fox_alphasynuclein <- aao_Fox %>% 
  filter(pathway ==  "alpha_synuclein") %>% 
  drop_na()

#Run meta-analysis using random-effects model
meta <- metagen(Coeff,
                se,
                data = aao_Fox_alphasynuclein,
                studlab = dataset,
                fixed = TRUE,
                random = TRUE,
                prediction = TRUE,
                sm = "SMD") #all effect sizes are in betas


