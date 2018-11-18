##################################proprocess annotation###############
##~~~~~~~~~~~~~~~~~~~~~~~~~~~gff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smugff <- read.delim('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.gff', skip = 3, header = FALSE, stringsAsFactor = FALSE)

names <- sapply(smugff[, 9], function(x) {
  eachNA <- unlist(strsplit(x, split = ';', fixed = TRUE))
  eachN <- eachNA[1]
  eachA <- eachNA[2]
  eachN <- unlist(strsplit(eachN, split = '=', fixed = TRUE))[2]
  eachA <- unlist(strsplit(eachA, split = '=', fixed = TRUE))[2]
  return(c(eachN, eachA))
})

names <- t(names)
colnames(names) <- c('geneNames', 'Anno')
smugff <- cbind(smugff, names)
smugff$loc <- paste(smugff$V4, smugff$V5, sep = '..')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~add ptt/rnt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smuptt <- read.delim('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.ptt', skip = 2, stringsAsFactor = FALSE)
smurnt <- read.delim('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.rnt', skip = 2, stringsAsFactor = FALSE)
smuannot <- rbind(smuptt, smurnt)

## deal with slash
slashLogic <- smuannot[, 'Gene'] == '-'
smuannot[slashLogic, 'Gene'] <- smuannot[slashLogic, 'Synonym']

smumerge <- merge(smugff, smuannot, by.x = 'loc', by.y = 'Location')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(smumerge, file = '/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/smumerge.RData')
######################################################################

###########################merge DEG list####################
library('magrittr')

setwd('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results')

load('smumerge.RData')

deg4h <- read.csv('4h_NC_004350_transcripts.txt', sep = '\t') %>%
  merge(smumerge, ., by.x = 'Synonym', by.y = 'Synonym')

deg24h <- read.csv('24h_NC_004350_transcripts.txt', sep = '\t') %>%
  merge(smumerge, ., by.x = 'Synonym', by.y = 'Synonym')

deg <- deg4h[, c(1, 26, 12, 6, 7, 13, 15, 19, 28:30, 36:38)]
deg <- cbind(deg, deg24h[, c(28:30, 36:38)])
## deg$log2FC <- log2((deg[, 'RPKM.2'] + 1)/(deg[, 'RPKM.1'] + 1))

colnames(deg) <- c('GeneID', 'Symbol', 'Names', 'Start', 'End', 'Product', 'Length', 'COG', paste0('RawCountControl4h', 1:3), paste0('RawCountMutant4h', 1:3), paste0('RawCountControl24h', 1:3), paste0('RawCountMutant24h', 1:3))

write.csv(deg, 'deg.csv')
#############################################################


##########################DEGseq2#################################
library('DESeq2')
library('ggplot2')
library('directlabels')
library('genefilter')
library('pheatmap')

setwd('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results')

##~~~~~~~~~~~~~~~~~~~~~~~~build target and DEGlist object~~~~~~~~~~~~~~
deg <- read.csv('deg.csv', row.names = 1, stringsAsFactor = FALSE)
## htCountSelect <- deg[, c(9:14)]
htCountSelect <- deg[, c(15:20)]
rownames(htCountSelect) <- deg[, 1]

targets <- data.frame(Group = factor(c('WT', 'WT', 'WT', 'deltasrtA', 'deltasrtA', 'deltasrtA')), Sample = paste0('smu', 1:6))
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
anno <- deg[match(rownames(res), deg[, 1]), 1:8]
res <- cbind(anno, res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
res <- res[order(res[, 'padj']), ]
write.csv(res, file = '/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/degseq24h_whole.csv')

## padj < 0.01 & |log2FC| > 1.5
sigLogic <- res$padj < 0.01 & abs(res$log2FoldChange) > log(1.5)
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
write.csv(resSig, file = '/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/degseq24h_DEG.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~
topNum <- nrow(resSig)
heatmapCount <- resSig[1:topNum, 9:14]
heatmapCount <- apply(heatmapCount, 1:2, as.numeric)
rownames(heatmapCount) <- resSig[1:topNum, 'Names']

annoCol <- data.frame(Group = colData(glioPR)[, 1])
row.names(annoCol) <- rownames(colData(glioPR))
annoColor <- list(Group = c(WT = '#00C19F', deltasrtA = '#F8766D'))
## annoRow = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(30, 30, 40))))
## rownames(annoRow) <- rownames(heatmapCount)
pdf('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/DESeq2_heatmap_whole.pdf')
pheatmap(heatmapCount, annotation_col = annoCol, annotation_colors = annoColor, fontsize=12, fontsize_row=2, annotation_legend = TRUE)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(glioPR, normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(glioPR,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
pdf('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/DESeq2_PCA.pdf')
groupCol <- c('#00C19F', '#F8766D')
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = groupCol) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
