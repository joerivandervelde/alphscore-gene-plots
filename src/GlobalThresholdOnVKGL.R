################################
# Install packages (only once) #
################################
#install.packages('ggplot2')
#install.packages("cutpointr")
#install.packages('toprdata')
#install.packages('seqminer')

#################
# Load packages #
#################
library(ggplot2)
library(cutpointr)
library(toprdata)
library(seqminer)

#######################################
# Set your working dir and file paths #
#######################################
setwd("/Users/joeri/git/alphscore-gene-plots/plots/")
# AlphScores, download from https://zenodo.org/record/6288139
alphScoreLoc <- "/Applications/AlphScore/AlphScore_final.tsv.gz"
# VKGL variant classifications, download from https://vkgl.molgeniscloud.org
vkglLoc <- "/Users/joeri/VKGL/VKGL-releases/VKGL_public_consensus_july2023.tsv"

####################################################
# Retrieve VKGL classifications and add AlphScores #
####################################################
vkgl <- read.table(file=vkglLoc, sep = '\t',header = TRUE)
uniqGenes <- unique(vkgl$gene)
results <- data.frame()
for(i in 1:length(uniqGenes)){
  geneName <- uniqGenes[i]
  cat(paste("Working on ",geneName," (",i," of ",length(uniqGenes),")\n", sep=""))
  geneCoords <- subset(ENSGENES, gene_symbol==geneName)
  if(nrow(geneCoords)==0){
    next
  }
  geneChr <- gsub("chr","", geneCoords$chrom)
  geneTabix <- paste(geneChr, paste(geneCoords$gene_start, geneCoords$gene_end, sep="-"), sep=":")
  geneAlphScoreData <- tabix.read(alphScoreLoc, geneTabix)
  alphScores <- read.table(text=geneAlphScoreData, sep = '\t',header = FALSE, fileEncoding = "UTF-16LE", col.names = c("chr","pos_1-based","ref","alt","aaref","aaalt","rs_dbSNP","hg19_chr","hg19_pos_1-based","ID","genename","Uniprot_acc_split","Uniprot_acc","HGVSp_VEP_split","HGVSp_VEP","CADD_raw","REVEL_score","DEOGEN2_score","b_factor","SOLVENT_ACCESSIBILITY_core","in_gnomad_train","in_clinvar_ds","AlphScore","glm_AlphCadd","glm_AlphRevel","glm_RevelCadd","glm_AlphRevelCadd","glm_AlphDeogen","glm_CaddDeogen","glm_DeogenRevel","glm_AlphDeogenRevel","glm_AlphCaddDeogen","glm_CaddDeogenRevel"))
  vkglGeneWithAlph <- merge(alphScores, vkgl, by.x = c("hg19_pos_1.based","ref","alt"), by.y = c("start","ref","alt"))
  vkglGeneWithAlph <- vkglGeneWithAlph[,c("AlphScore","classification")] # Leave out "pos_1.based","ref","alt", etc
  vkglGeneWithAlph <- subset(vkglGeneWithAlph, classification != "VUS")
  #cat(paste("Adding ",nrow(vkglGeneWithAlph)," variants\n", sep=""))
  results <- rbind(results, vkglGeneWithAlph)
}

####################################################
# Determine optimal threshold using Youden's Index #
####################################################
opt_cut <- cutpointr(results, AlphScore, classification, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(results[results$classification=="LP",'AlphScore'] >= youdenIndex)
fp <- sum(results[results$classification=="LB",'AlphScore'] >= youdenIndex)
tn <- sum(results[results$classification=="LB",'AlphScore'] < youdenIndex)
fn <- sum(results[results$classification=="LP",'AlphScore'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- opt_cut$sensitivity*100
spec <- opt_cut$specificity*100
cat(paste("The optimal AlphScore threshold based on VKGL variant classifications is ",round(youdenIndex,5)," with PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))

# VKGL July 2023 release:
# The optimal AlphScore threshold based on VKGL variant classifications is 0.79734 with PPV 38%, NPV 90%, sensitivity 66% and specificity 75%.

write.table(results, sep="\t",file="alphscore-vkgl-results.txt", quote=FALSE, row.names =FALSE)
