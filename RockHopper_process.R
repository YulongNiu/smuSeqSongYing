##################################proprocess annotation###############
##~~~~~~~~~~~~~~~~~~~~~~~~~~~gff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smugff <- read.delim('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.gff', skip = 3, header = FALSE, stringsAsFactor = FALSE)

names <- sapply(smugff[, 9], function(x) {
  eachNA <- unlist(strsplit(x, split = ';', fixed = TRUE))
  eachN <- eachNA[1]
  eachA <- eachNA[2]
  eachN <- unlist(strsplit(eachN, split = '=', fixed = TRUComputational analysis of bacterial RNA-seq dataE))[2]
  eachA <- unlist(strsplit(eachA, split = '=', fixed = TRUE))[2]
  return(c(eachN, eachA))
})

names <- t(names)
colnames(names) <- c('geneNames', 'Anno')
smugff <- cbind(smugff, names)
smugff$loc <- paste(smugff$V4, smugff$V5, sep = '..')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~add ptt/rnt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smuptt <- read.delim('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.ptt', skip = 2, stringsAsFactor = FALSE)
smurnt <- read.delim('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.rnt', skip = 2, stringsAsFactor = FALSE)
smuannot <- rbind(smuptt, smurnt)

## deal with slash
slashLogic <- smuannot[, 'Gene'] == '-'
smuannot[slashLogic, 'Gene'] <- smuannot[slashLogic, 'Synonym']

smumerge <- merge(smugff, smuannot, by.x = 'loc', by.y = 'Location')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(smumerge, file = '/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/smumerge.RData')
######################################################################

############################map wig############################
MapWig <- function(wig, start, end) {
  ## INPUT: 'wig' is a numeric vector. 'start' and 'end' indicate the start and end position.
  ## OUTPUT: the count number

  return(sum(wig[start:end]))
}

read.wig <- function(wigpath) {
  wig <- read.table(wigpath, skip = 2, stringsAsFactors = FALSE, header = FALSE)
  return(wig[, 1])
}

setwd('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/genomeBrowserFiles/')

library('magrittr')
library('foreach')
library('doMC')

filename <- dir(pattern = 'v1')
tmp1 <- foreach (i = 1:length(filename), .combine = cbind) %do% {
  return(read.wig(filename[i]))
}
###############################################################


###########################merge DEG list####################
setwd('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results')

load('smumerge.RData')
deg <- read.csv('NC_004350_transcripts.txt', sep = '\t')
deg <- merge(smumerge, deg, by.x = 'Synonym', by.y = 'Synonym')

deg$log2FC <- log2((deg[, 'RPKM.2'] + 1)/(deg[, 'RPKM.1'] + 1))
deg <- deg[, c(1, 26, 6, 7, 13, 15, 19, 28:30, 34, 36:38, 42, 46, 44, 45)]
colnames(deg) <- c('GeneID', 'GeneName', 'Start', 'End', 'Product', 'Length', 'COG',paste0('RawCountControl', 1:3), 'RPKMControl', paste0('RawCountMutant', 1:3), 'RPKMMutant', 'log2FC', 'p-value', 'q-value')

write.csv(deg, 'deg.csv')
#############################################################


##########################DEGseq2#################################
library('DESeq2')
library('ggplot2')
library('directlabels')
library('genefilter')
library('PoiClaClu')
library('pheatmap')

setwd('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results')


##~~~~~~~~~~~~~~~~~~~~~~~~build target and DEGlist object~~~~~~~~~~~~~~
deg <- read.csv('deg.csv', row.names = 1, stringsAsFactor = FALSE)
htCountSelect <- deg[, c(8:10, 12:14)]
rownames(htCountSelect) <- deg[, 1]

targets <- data.frame(Group = factor(c('WT', 'WT', 'WT', 'Mutant', 'Mutant', 'Mutant')), Sample = paste0('smu', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(htCountSelect) <- rownames(targets)
glioPR <- DESeqDataSetFromMatrix(countData = htCountSelect, colData = targets, design = ~Group)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEG analysis~~~~~~~~~~~~~~~~~~~~~~~~
glioPR <- glioPR[rowSums(counts(glioPR)) > 1, ]
glioPR <- DESeq(glioPR)
## count transformation
rld <- rlog(glioPR)
vst <- varianceStabilizingTransformation(glioPR)
resRaw <- results(glioPR)
resRaw[, 2] <- -resRaw[, 2]
summary(resRaw)
res <- cbind(as.matrix(mcols(glioPR)[, 1:10]), assay(rld))
anno <- deg[match(rownames(res), deg[, 1]), 1:7]
res <- cbind(anno, res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
write.csv(res, file = '/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/degseq2.csv')
resSig <- res[order(as.numeric(res[, 'qval'])), ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
