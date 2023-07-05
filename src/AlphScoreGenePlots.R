################################
# Install packages (only once) #
################################
#install.packages('seqminer')
#install.packages('toprdata')
#install.packages('ggplot2')
#install.packages('stringr')
#install.packages('scales')
#install.packages("cutpointr")
#install.packages('plyr')
#install.packages('dplyr')

#################
# Load packages #
#################
library(seqminer)
library(toprdata)
library(ggplot2)
library(stringr)
library(scales)
library(cutpointr)
library(plyr)
library(dplyr)

##############################
# Set your gene name to plot #
##############################
geneName <- "COL3A1"

#######################################
# Set your working dir and file paths #
#######################################
setwd("/Users/joeri/git/alphscore-gene-plots/plots/")
# AlphScores, download from https://zenodo.org/record/6288139
alphScoreLoc <- "/Applications/AlphScore/AlphScore_final.tsv.gz"
# VKGL variant classifications, download from https://vkgl.molgeniscloud.org
vkglLoc <- "/Users/joeri/VKGL-releases/VKGL_public_consensus_apr2023.tsv"
# ClinVar variant classifications, download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
clinVarLoc <- "/Applications/ClinVar/clinvar_20230702.vcf.gz"

######################################
# Retrieve gene and exon coordinates #
######################################
geneCoords <- subset(ENSGENES, gene_symbol==geneName)
geneChr <- gsub("chr","", geneCoords$chrom)
exonCoords <- subset(ENSEXONS, gene_symbol==geneName)
exonStartCoords <- str_split(exonCoords$exon_chromstart, ",")
exonEndCoords <- str_split(exonCoords$exon_chromend, ",")
exons <- data.frame(exonStart = as.numeric(unlist(exonStartCoords)), exonEnd = as.numeric(unlist(exonEndCoords)))

#######################
# Retrieve AlphScores #
#######################
geneTabix <- paste(geneChr, paste(geneCoords$gene_start, geneCoords$gene_end, sep="-"), sep=":")
geneAlphScoreData <- tabix.read(alphScoreLoc, geneTabix)
alphScores <- read.table(text=geneAlphScoreData, sep = '\t',header = FALSE, fileEncoding = "UTF-16LE", col.names = c("chr","pos_1-based","ref","alt","aaref","aaalt","rs_dbSNP","hg19_chr","hg19_pos_1-based","ID","genename","Uniprot_acc_split","Uniprot_acc","HGVSp_VEP_split","HGVSp_VEP","CADD_raw","REVEL_score","DEOGEN2_score","b_factor","SOLVENT_ACCESSIBILITY_core","in_gnomad_train","in_clinvar_ds","AlphScore","glm_AlphCadd","glm_AlphRevel","glm_RevelCadd","glm_AlphRevelCadd","glm_AlphDeogen","glm_CaddDeogen","glm_DeogenRevel","glm_AlphDeogenRevel","glm_AlphCaddDeogen","glm_CaddDeogenRevel"))

####################################
# Retrieve ClinVar classifications #
####################################
geneClinVarData <- tabix.read(clinVarLoc, geneTabix)
geneClinVarDF <- read.table(text=geneClinVarData, sep = '\t',header = FALSE, fileEncoding = "UTF-16LE", col.names = c("chr","pos","id","ref","alt","qual","filter","info"))
# Remove any conflicts and records without a classification
clinvar <- subset(geneClinVarDF, !grepl("Conflicting_interpretations_of_pathogenicity", geneClinVarDF$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_assertion_provided", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("drug_response", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_assertion_criteria_provided", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_interpretation_for_the_single_variant", clinvar$info, fixed=TRUE))
# Relabel classifications into 3 groups
clinvar$classification[grepl("benign", clinvar$info, ignore.case = TRUE)] <- "LB/B"
clinvar$classification[grepl("uncertain_significance", clinvar$info, ignore.case = TRUE)] <- "VUS"
clinvar$classification[grepl("pathogenic", clinvar$info, ignore.case = TRUE)] <- "LP/P"
clinvar$Source <- "ClinVar"
# there should be 0 records left without a classification, this is a way to check:
# subset(clinvar, clinvar$classification==FALSE)

#################################
# Retrieve VKGL classifications #
#################################
vkglDataFrame <- read.table(file=vkglLoc, sep = '\t',header = TRUE)
vkgl <- subset(vkglDataFrame, gene==geneName)
vkgl$classification <- revalue(vkgl$classification, c("VUS"="VUS", "LB"="LB/B", "LP"="LP/P"))
vkgl$Source <- "VKGL"

##########################################
# Merge ClinVar and VKGL with AlphScores #
##########################################
clinvarWithAlph <- merge(alphScores, clinvar, by.x = c("pos_1.based","ref","alt"), by.y = c("pos","ref","alt"))
clinvarWithAlph <- clinvarWithAlph[,c("pos_1.based","ref","alt","AlphScore","classification","Source")]
# VKGL is build 37 i.e. hg19 data, luckily we can still merge using the hg19_pos_1.based column
vkglWithAlph <- merge(alphScores, vkgl, by.x = c("hg19_pos_1.based","ref","alt"), by.y = c("start","ref","alt"))
vkglWithAlph <- vkglWithAlph[,c("pos_1.based","ref","alt","AlphScore","classification","Source")]
# Combine all together and order by classification label for a better plot (LP/P on top, then LB/B, then VUS)
variants <- rbind(vkglWithAlph,clinvarWithAlph)
variants <- variants %>% arrange(factor(classification, levels = c("VUS","LB/B","LP/P")))

####################################################
# Determine optimal threshold using Youden's Index #
####################################################
cutpointDF <- subset(variants, classification != "VUS")
opt_cut <- cutpointr(cutpointDF, AlphScore, classification, direction = ">=", pos_class = "LP/P", neg_class = "LB/B", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(cutpointDF[cutpointDF$classification=="LP/P",'AlphScore'] >= youdenIndex)
fp <- sum(cutpointDF[cutpointDF$classification=="LB/B",'AlphScore'] >= youdenIndex)
ppv <- 100 *tp/(tp+fp)
sens <- opt_cut$sensitivity*100

#########################
# Determine plot window #
#########################
xmin <- min(exons$exonStart)
xmax <- max(exons$exonEnd)
ymin <- min(alphScores$AlphScore)
ymax <- max(alphScores$AlphScore)

########################
# Create and save plot #
########################
ggplot() +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
  geom_rect(data = exons, aes(xmin = exonStart, xmax = exonEnd, ymin = ymin, ymax = ymax), linetype = 0, fill="lightgray", alpha = 1) +
  geom_point(data = alphScores, aes(x=pos_1.based, y=AlphScore), alpha=1.0, size = 0.5, colour="black") +
  geom_point(data = variants, aes(x=pos_1.based, y=AlphScore, colour=classification, shape=Source), alpha=1.0, size = 2, stroke = 1) +
  geom_hline(yintercept = youdenIndex) +
  scale_colour_manual(name = "Classification", values = c("LB/B" = "green","VUS" = "darkgray","LP/P" = "red")) +
  scale_shape_manual(name = "Source", values = c("ClinVar" = 6,"VKGL" = 2)) +
  scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
  xlab(paste("",geneName," at GRCh38 chr",geneChr,":", xmin, "-",xmax,". Lightgray: exons, black dots: scored potential missense variants.", sep="")) +
  ylab("AlphScore (Schmidt et al., Bioinformatics, May 2023)") +
  ggtitle(paste("AlphScores for ",geneName,". At a threshold of ",round(youdenIndex, 2), " the PPV is ",round(ppv),"% and the sensitivity is ",round(sens),"%.",sep=""))
ggsave(paste(geneName,".png",sep=""), width=9, height=5)
