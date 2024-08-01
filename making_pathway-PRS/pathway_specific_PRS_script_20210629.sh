### PATHWAY-SPECIFIC PRS ###
#Created: 03/05/2021
#Last updated: 03/05/2021
#Created by: Manuela Tan
#WD: /Users/manuela/Documents/Work/pathway_specific_PRS/

### Pathway of interest ###
	
	#Alpha-synuclein pathway (Bandres-Ciga et al., 2020)
	#Adaptive immunity (Bandres-Ciga et al., 2020; Holmans et al., 2013)
	#Innate immunity (Bandres-Ciga et al., 2020; Holmans et al., 2013)
	#Lysosome pathway (Bandres-Ciga et al., 2020; Robak et al., 2017)
	#Endocytic membrane-trafficking pathway (Bandres-Ciga et al., 2019)
	#Mitochondrial pathway (Billingsley et al., 2019)
	#Microglial open chromatin SNPs (Andersen et al., 2021)
	#Monocyte open chromatin SNPs (Andersen et al., 2021)

### Get pathways of interest from Molecular Signature Database (c2.cp.v7.0.symbols.gmt) ###

	#Alpha-synuclein pathway
	grep SYNUCLEIN c2.cp.v7.0.symbols.gmt > alpha_synuclein_genes.txt

	#Adaptive immunity
	grep REACTOME_ADAPTIVE_IMMUNE_SYSTEM c2.cp.v7.0.symbols.gmt > adaptive_immune_genes.txt

	#Innate immunity
	grep REACTOME_INNATE_IMMUNE_SYSTEM c2.cp.v7.0.symbols.gmt > innate_immune_genes.txt

	#There is also a pathway called REACTOME_REGULATION_OF_INNATE_IMMUNE_RESPONSES_TO_CYTOSOLIC_DNA
	#But all these genes are captured in the REACTOME_INNATE_IMMUNE_SYSTEM gene list

	#Lysosome
	grep KEGG_LYSOSOME c2.cp.v7.0.symbols.gmt > lysosome_genes.txt


### Get gene bp positions from BioMart GRCh37 ###
	#http://grch37.ensembl.org/biomart/martview/
	#Database = Ensembl Genes 103
	#Dataset = Human genes (GRCh37.p13)
	#Filters > Gene 
		#Input external references ID list
		#Gene name
	#Attributes > Gene
		#Gene stable ID
		#Gene name
		#Chromosome/scaffold name
		#Gene start (bp)
		#Gene end (bp)
	#Then go to Results
	#Export as CSV file - keep unique results only
	#I named the file $PATHWAY_genes_positions.txt

	#Read into R and format for plink


### For microglia and monocyte open chromatin regions, use files from Corces and Nott papers ###
	
	#files in /Users/manuela/Documents/Work/pathway_specific_PRS/pathway_snplists/Corces_Nott_ATACseq
	#For microglia from Nott: PU1_optimal_peak_IDR_ENCODE.ATAC.bed
	#For monocytes from Corces: atac_monocytes_Corces_c50.bed
	#Need to reformat these files slightly to extract in plink format
	#See R script

### In plink, extract SNPs in each pathway list ###

### Liftover to hg38 for AMP-PD data ###
	#From hg19 to hg38
	
	cd ./liftover_hg38 

	#Download chain file
	rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz .

	gunzip hg19ToHg38.over.chain.gz 


	#First do some formatting in R - in pathway_specific_PRS_script_20210629.R
	#Run liftOver - only on ATACseq microglia and monocyte pathways
	for PATHWAY in microglia monocytes
	do
	/Users/manuela/Documents/Work/software/liftOver/liftOver -multiple ../outputs/hg38_liftover/"$PATHWAY".txt hg19ToHg38.over.chain output_"$PATHWAY".bed unlifted_"$PATHWAY".bed
	done

	#For the pathways that are gene lists - will use BioMart to get hg38 coordinates (follow BioMart steps above just using hg38)
	#Only use lifted over coordinates for the pathways that are not gene lists - microglial open chromatin and monocyte open chromatin




