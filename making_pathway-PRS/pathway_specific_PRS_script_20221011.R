##### PATHWAY SPECIFIC PRS #####

#---Load packages---####
library(tidyverse)
library(data.table)
library(readxl)

#---Read in pathways from MsigDB data, reformat and export---####

#Make list of pathways
#These are just the ones that we do not have chr:bp positions for
#The other pathways are in the correct format for extracting from plink
pathway_list <- c("adaptive_immune", "alpha_synuclein", "innate_immune", "lysosome")


#There are some genes where the chromosome/scaffold is not a chr number but a weird patch
#https://www.biostars.org/p/106355/
#Need to map these to actual chromosomes
#Patch data downloaded from 
#ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/PATCHES/alt_scaffolds/alt_scaffold_placement.txt
#File called alt_scaffold_placement.txt

alt_scaffold <- fread("./pathway_snplists/alt_scaffold_placement.txt")

#Select just relevant columns
alt_scaffold <- alt_scaffold %>% 
  select(alt_scaf_name, parent_name, parent_start, parent_stop, region_name)

for (nm in 1:length(pathway_list)) {
  
  #Read in data
  data <- read.csv(paste("pathway_snplists/", pathway_list[nm], "_genes_positions.txt", sep = ""))

  #Remove duplicates - keep the one with the numeric Chromosome scaffold name
  data_unique <- data %>% 
    arrange(Chromosome.scaffold.name) %>% 
    distinct(Gene.name, .keep_all = TRUE)
  
  #Map patch names
  data_unique_merged <- data_unique %>% 
    left_join(alt_scaffold, by = c("Chromosome.scaffold.name" = "alt_scaf_name")) 
  
  #Separate Chromosome.scaffold.name to extract HSCHR section
  data_unique_merged <- data_unique_merged %>%
    mutate(scaffold_name = Chromosome.scaffold.name) %>% 
    separate(scaffold_name, into = c("HSCHR", NA, NA), sep = "_")
  #Ignore warnings
  
  #If scaffold name contains "CHR", then make new column called chr which removes all letters
  data_unique_merged <- data_unique_merged %>% 
    mutate(chr_numeric = as.numeric(as.character(Chromosome.scaffold.name))) %>%
    mutate(chr_hschr = ifelse(grepl("CHR", HSCHR, fixed = TRUE) == TRUE, 
                              (gsub("[^0-9.-]", "", HSCHR)), 
                              "no")) %>%
    mutate(chr_final = ifelse(!is.na(chr_numeric), chr_numeric,
                              ifelse(is.na(chr_numeric) & (Chromosome.scaffold.name == "X" | Chromosome.scaffold.name == "Y"), Chromosome.scaffold.name,
                                     ifelse(is.na(chr_numeric) & chr_hschr!="no", chr_hschr,
                                            ifelse(is.na(chr_numeric) & chr_hschr == "no", parent_name, NA)))))
  
  
  #Check any missing final chromosome
  print(pathway_list[nm])
  
  print(data_unique_merged %>% 
    filter(is.na(chr_final))) %>% 
    summarise(count = n())
  
  #Export for plink format
  export <- data_unique_merged %>% 
    select(chr_final, Gene.start..bp., Gene.end..bp., Gene.name)
  
  #Filename with pathway name
  filename <- paste("./outputs/", pathway_list[nm], "_extract.txt", sep = "")
  
  #Export in tab separated file
  write.table(export, filename, 
              sep = "\t", quote = F, col.names = F, row.names = F)
}




#---Read in microglia ATACseq data, reformat and export---####

#Read in microglia ATACseq data from Nott
microglia <- read.table("pathway_snplists/Corces_Nott_ATACseq/PU1_optimal_peak_IDR_ENCODE.ATAC.bed")

#Make new column for chromosome without "chr" text
microglia <- microglia %>% 
  mutate(chr = gsub("chr", "", V1))

microglia_export <- microglia %>% 
  select(chr, V2, V3, V1)

#Export as tab separated file
write.table(microglia_export, "./outputs/microglia_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)


#---Read in monocyte ATACseq data, reformat and export---####

#Read in monocyte ATACdeq data from Corces
monocyte <- fread("pathway_snplists/Corces_Nott_ATACseq/atac_monocytes_Corces_c50.bed")

#Select just chr, start and end bp columns
monocyte_selected <- monocyte[,c(1:3)]

#Rename column 1
colnames(monocyte_selected) <- c("chr", "startbp", "endbp")

#Make new column for chromosome
monocyte_selected <- monocyte_selected %>% 
  mutate(chr_new = gsub("chr", "", chr))

#Format for export
monocyte_export <- monocyte_selected %>% 
  select(chr_new, startbp, endbp, chr)

#Export as tab separated file
write.table(monocyte_export, "./outputs/monocytes_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)

#---Read in mitochondria gene lists from Billingsley, reformat and export---####

#Read in gene lists from Excel spreadsheet (data from Supplementary Tables 3 and 4)

mitochondria_primary <- read_xlsx("./pathway_snplists/mitochondria_gene_lists_Billingsley2019.xlsx",
                                  sheet = 1)

mitochondria_secondary <- read_xlsx("./pathway_snplists/mitochondria_gene_lists_Billingsley2019.xlsx",
                                    sheet = 2)

#These are already in the right format for plink extract

#Export as tab separated file
write.table(mitochondria_primary, "./outputs/mitochondria_primary_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)

write.table(mitochondria_secondary, "./outputs/mitochondria_secondary_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)


#Check overlap between gene lists
overlap_mitochondria <- mitochondria_primary %>% 
  inner_join(mitochondria_secondary, by = "Gene")

#Make a merged mitochondria list
mitochondria_all <- rbind(mitochondria_primary, mitochondria_secondary)

#Remove duplicates
mitochondria_all_unique <- mitochondria_all %>% 
  distinct(Gene, .keep_all = TRUE)

#Export as tab separated file
write.table(mitochondria_all_unique, "./outputs/mitochondria_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)

#Write just the list of mitochondria genes - so we can get hg38 positions from BioMart
mitochondria_genes <- mitochondria_all_unique %>% 
  select(Gene)

write.table(mitochondria_genes, "./pathway_snplists/mitochondria_genes.txt",
            sep = "\t", quote = F, col.names = F, row.names = F)


#---Read in endoocytosis gene lists sent by Sara---####

#Read in gene list
endocytosis <- read.table("pathway_snplists/endocytosis.txt")

#This is already in the right format for plink extract

#Export as tab separated file
write.table(endocytosis, "./outputs/endocytosis_extract.txt", 
            sep = "\t", quote = F, col.names = F, row.names = F)

#Write just the list of endocytosis genes - so we can get hg38 positions from BioMart
endocytosis_genes <- endocytosis %>% 
  select(V4)

write.table(endocytosis_genes, "./pathway_snplists/endocytosis_genes.txt",
            sep = "\t", quote = F, col.names = F, row.names = F)



#---Read in final files and format for liftover---####

#Just microglia and monocytes ATACseq
pathways_liftover <- c("microglia", "monocytes")

for (nm in 1:length(pathways_liftover)) {
  
  pathway <- pathways_liftover[nm]
  
  #Read in formatted pathway extract lists in hg19
  data <- read.table(paste("outputs/", pathway, "_extract.txt", sep = ""), header = F)
  
  #Label column names
  colnames(data) <- c("chr", "start_bp", "end_bp", "gene")
  
  #For liftover, regions should be in chrN:start-end formats
  data <- data %>% 
    mutate(liftover_chr = paste("chr", chr, sep = ""))
  
  #Export columns for liftover
  data_export <- data %>% 
    select(liftover_chr, start_bp, end_bp, gene)
  
  filename <- paste("./outputs/hg38_liftover/", pathway, ".txt", sep = "")
  
  write.table(data_export, filename, quote = F, row.names = F, col.names = F, sep = "\t")
  
}

#Run liftover in terminal

#Read in liftover files
for (nm in 1:length(pathways_liftover)) {
  
  pathway <- pathways_liftover[nm]
  
  output <- read.table(paste("./liftover_hg38/output_", pathway, ".bed", sep = ""))
  
  #Format for export
  output <- output %>% 
    mutate(chr = gsub("chr", "", V1))
  
  hg38_export <- output %>% 
    select(chr, V2, V3)
  
  filename <- paste("./outputs/hg38_liftover/", pathway, "_extract_hg38.txt", sep = "")
  
  write.table(hg38_export, filename, quote = F, col.names = F, row.names = F)
  
}
  
#---Read in hg38 data, reformat and export---####

#Using BioMart data
#For all pathways except for the microglia and monocyte ATACseq
pathways_hg38 <- c("adaptive_immune", "alpha_synuclein", "innate_immune", "lysosome", "mitochondria", "endocytosis")

#There are some genes where the chromosome/scaffold is not a chr number but a weird patch
#https://www.biostars.org/p/106355/
#Need to map these to actual chromosomes
#Patch data downloaded from 
#https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38.p3/PATCHES/alt_scaffolds/alt_scaffold_placement.txt
#File called alt_scaffold_placement.txt

alt_scaffold <- fread("./pathway_snplists/alt_scaffold_placement_hg38.txt")

#Select just relevant columns
alt_scaffold <- alt_scaffold %>% 
  select(alt_scaf_name, parent_name, parent_start, parent_stop, region_name)

for (nm in 1:length(pathways_hg38)) {
  
  #Read in data
  data <- read.csv(paste("pathway_snplists/", pathways_hg38[nm], "_genes_positions_hg38.txt", sep = ""))
  
  #Remove duplicates - keep the one with the numeric Chromosome scaffold name
  #These bp coordinates match gnomAD positions
  data_unique <- data %>% 
    arrange(Chromosome.scaffold.name) %>% 
    distinct(Gene.name, .keep_all = TRUE)
  
  #Map patch names
  data_unique_merged <- data_unique %>% 
    left_join(alt_scaffold, by = c("Chromosome.scaffold.name" = "alt_scaf_name")) 
  
  #Separate Chromosome.scaffold.name to extract HSCHR section
  data_unique_merged <- data_unique_merged %>%
    mutate(scaffold_name = Chromosome.scaffold.name) %>% 
    separate(scaffold_name, into = c("CHR", "HSCHR", NA, NA, NA), sep = "_")
  #Ignore warnings
  
  #If scaffold name contains "CHR", then make new column called chr which removes all letters
  data_unique_merged <- data_unique_merged %>% 
    mutate(chr_numeric = as.numeric(as.character(Chromosome.scaffold.name))) %>%
    mutate(chr_hschr = ifelse(grepl("CHR", HSCHR, fixed = TRUE) == TRUE, 
                              (gsub("[^0-9.-]", "", HSCHR)), 
                              "no")) %>%
    mutate(chr_final = ifelse(!is.na(chr_numeric), chr_numeric,
                              ifelse(is.na(chr_numeric) & (Chromosome.scaffold.name == "X" | Chromosome.scaffold.name == "Y"), Chromosome.scaffold.name,
                                     ifelse(is.na(chr_numeric) & chr_hschr!="no", chr_hschr,
                                            ifelse(is.na(chr_numeric) & chr_hschr == "no", parent_name, NA)))))
  
  
  #Check any missing final chromosome
  print(pathways_hg38[nm])
  
  print(data_unique_merged %>% 
          filter(is.na(chr_final))) %>% 
    summarise(count = n())
  
  #Export for plink format
  export <- data_unique_merged %>% 
    select(chr_final, Gene.start..bp., Gene.end..bp., Gene.name)
  
  #Filename with pathway name
  filename <- paste("./outputs/hg38_gene_coords/", pathways_hg38[nm], "_extract_hg38.txt", sep = "")
  
  #Export in tab separated file
  write.table(export, filename, 
              sep = "\t", quote = F, col.names = F, row.names = F)
}


#---Look at overlap between gene lists---####


pathways <- c("adaptive_immune", "alpha_synuclein", "innate_immune", "lysosome", "mitochondria", "endocytosis", "microglia", "monocytes")

for (i in 1:length(pathways)) {
  
  #Read in data
  df <- read.table(paste("./outputs/", pathways[i], "_extract.txt", sep = ""))
  colnames(df) <- c("chr", "bp_start", "bp_end", "gene")
  
  #Add a column for pathway
  df <- df %>% 
    mutate(pathway = pathways[i])
  
  #Save as dataframe
  assign(pathways[i], df)
  
  rm(df)
  
}


#Make function to get the overlapping genes between pathway lists
#df is the dataframe to compare to all the others
#Should be in character format
#df="adaptive_immune"
get_overlap <- function(df){
  
  #Get list of other pathways
  other_pathways <- setdiff(pathways, df)
  
  overlapping_list <- list()
  
  for (i in 1:length(other_pathways)) {
    
    #Get first dataset for comparison
    compare_df <- get(other_pathways[i])
    
    #Get original dataset
    data <- get(df)
    
    #Get overlapping genes between original df and compare df
    overlapping <- data %>% 
      select(pathway, gene) %>%
      inner_join(compare_df, by = "gene")
    
    #Add to main dataframe
    overlapping_list[[i]] <- overlapping
    
  }
  
  overlapping_final <- do.call(rbind, overlapping_list)
  
}


#Look at overlapping genes between pathway lists (by gene)
adaptive_immune_overlap <- get_overlap("adaptive_immune")
write.table(adaptive_immune_overlap, "./outputs/adaptive_immune_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(adaptive_immune_overlap)/nrow(adaptive_immune)

alpha_synuclein_overlap <- get_overlap("alpha_synuclein")
write.table(alpha_synuclein_overlap, "./outputs/alpha_synuclein_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(alpha_synuclein_overlap)/nrow(alpha_synuclein)

endocytosis_overlap <- get_overlap("endocytosis")
write.table(endocytosis_overlap, "./outputs/endocytosis_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(endocytosis_overlap)/nrow(endocytosis)
endocytosis_overlap$gene[duplicated(endocytosis_overlap$gene)]

innate_immune_overlap <- get_overlap("innate_immune")
write.table(innate_immune_overlap, "./outputs/innate_immune_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(innate_immune_overlap)/nrow(innate_immune)

lysosome_overlap <- get_overlap("lysosome")
write.table(lysosome_overlap, "./outputs/lysosome_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(lysosome_overlap)/nrow(lysosome)

mitochondria_overlap <- get_overlap("mitochondria")
write.table(mitochondria_overlap, "./outputs/mitochondria_overlap.txt", quote = F, col.names = T,
            row.names = F, sep = "\t")

nrow(mitochondria_overlap)/nrow(mitochondria)


#Look at overlapping regions for microglia and monocyte open chromatin regions

  

#---Reformat gene lists for LDSC---####


pathways <- c("adaptive_immune", "alpha_synuclein", "innate_immune", "lysosome", "mitochondria", "endocytosis", "microglia", "monocytes")

for (i in 1:length(pathways)) {
  
  #Read in data
  df <- read.table(paste("./outputs/", pathways[i], "_extract.txt", sep = ""))
  colnames(df) <- c("chr", "bp_start", "bp_end", "gene")
  
  #Remove gene column
  df_selected <- df %>% 
    mutate(chr_new = paste("chr", chr, sep = "")) %>% 
    select(chr_new, bp_start, bp_end)
  
  filename <- paste("./outputs/for_LDSC/", pathways[i], ".bed", sep = "")
  
  #Export into output folder
  write.table(df_selected, filename,
              quote = F, col.names = F, row.names = F, sep = "\t")
}
