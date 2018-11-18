############################kallisto annotation##############
library('Biostrings')
library('stringr')
library('magrittr')

deg <- read.csv('/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/deg.csv', row.names = 1, stringsAsFactor = FALSE) %>%
  `[`(., , 1:8)

gff <- read.table('/home/Yulong/Biotools/RefData/smu/NC_004350.gff', skip = 3, header = FALSE, sep = '\t', stringsAsFactors = FALSE) %>%
  `[`(., , 9) %>%
  strsplit(., split = ';', fixed = TRUE) %>%
  sapply(., `[`, 1) %>%
  strsplit(., split = '=', fixed = TRUE) %>%
  sapply(., `[`, 2)

deg <- deg[order(deg$Names), ][rank(gff), ]

smucdna <- readBStringSet('/home/Yulong/Biotools/RefData/smu/NC_004350_cdna.fa')
names(smucdna) <- deg$GeneID

writeXStringSet(smucdna, '/home/Yulong/Biotools/RefData/smu/NC_004350_cdna_name.fa')
#############################################################

###########################k res####################
##~~~~~~~~~~~~~~~~~~~~~~~~~import~~~~~~~~~~~~~~~~~~~~~
library('tximport')
library('rhdf5')
library('DESeq2')
library('edgeR')
library('magrittr')

wd <- '/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results'

fmember = c('4h_sm_1', '4h_sm_2', '4h_sm_3', '4h_srta_1', '4h_srta_2', '4h_srta_3', '24h_sm_1', '24h_sm_2', '24h_sm_3', '24h_srta_1', '24h_srta_2', '24h_srta_3')

files <- file.path(wd, fmember, 'abundance.h5')
names(files) <- fmember
txi.kallisto <- tximport(files, type = 'kallisto', txOut = TRUE)

htCountSelect <- txi.kallisto
htCountSelect$abundance %<>% `[`(., , 1:6)
htCountSelect$counts %<>% `[`(., , 1:6)
htCountSelect$length %<>% `[`(., , 1:6)

targets <- data.frame(Group = factor(c('WT', 'WT', 'WT', 'deltasrtA', 'deltasrtA', 'deltasrtA')), Sample = paste0('smu', 1:6))
rownames(targets) <- paste0(targets$Group, c(1:3, 1:3))
colnames(htCountSelect$counts) <- rownames(targets)
glioPR <- DESeqDataSetFromTximport(htCountSelect, colData = targets, design = ~Group)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DESeq2 analysis~~~~~~~~~~~~~~~~~~~~
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

res <- cbind(res[, 11:16], data.frame(resRaw[, c(5, 6, 2)]))
res <- res[order(res[, 'padj']), ]
write.csv(res, file = '/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/degseq24h_whole.csv')

## padj < 0.01 & |log2FC| > 1
sigLogic <- res$padj < 0.01 & abs(res$log2FoldChange) > 1
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
write.csv(resSig, file = '/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/degseq24h_DEG.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~edgeR~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cts <- htCountSelect$counts
normMat <- htCountSelect$length
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))

## filtering
keep <- filterByExpr(y)
y <- y[keep, ]

group <- factor(rep(1:2, each = 3))
design <- model.matrix(~group)
y <- estimateDisp(y, design)

## quasi-likelihood F-tests
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit,coef=2)
res <- topTags(qlf, n = 1949)
res <- res$table
sigLogic <- res$FDR < 0.01 & abs(res$logFC) > 1
sigLogic[is.na(sigLogic)] <- FALSE
resSig <- res[sigLogic, ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################
