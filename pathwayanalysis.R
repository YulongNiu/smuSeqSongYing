#####################KEGG####################################
setwd('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results')

library('KEGGAPI')

smuPathRaw <- getKEGGPathGenes('smu')
smuPathRaw <- sapply(smuPathRaw, function(x) {
  eachID <- sapply(strsplit(x, split = ':', fixed = TRUE), '[[', 2)
  return(eachID)
})

smuIDs <- res[, 1]
smuKEGG <- lapply(smuPathRaw, function(x) {
  return(x[x %in% smuIDs])
})

save(smuKEGG, file = 'smuKEGG.RData')
#############################################################

##########################GO#################################
setwd('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results')

load('kallisto_results/smuGO.RData')

res <- read.csv('degseq4h_whole.csv', stringsAsFactor = FALSE)

smuIDs <- res[, 1]
smuGO <- lapply(smuGO, function(x) {
  return(x[x %in% smuIDs])
})

save(smuGO, file = 'kallisto_results/smuGO.RData')
#############################################################

##########################GO analysis###########################
setwd('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('magrittr')

registerDoMC(8)

load('smuGO.RData')
load('smuKEGG.RData')
res <- read.csv('degseq4h_whole.csv', stringsAsFactor = FALSE)

## remove 0 terms
smuGO %<>% `[`(sapply(smuGO, length) > 0)
smuKEGG %<>% `[`(sapply(smuKEGG, length) > 0)

## padj < 0.01 & |log2FC| > 1
degVecLogic <- res$padj < 0.01 & abs(res$log2FoldChange) > log2(2)
degVecLogic[is.na(degVecLogic)] <- FALSE
degVec <- as.integer(degVecLogic)
names(degVec) <- res$GeneID

pwf <- nullp(degVec, bias.data = res$Length)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GOMat <- foreach(i = 1:length(smuGO), .combine = rbind) %dopar% {
  eachMat <- cbind(smuGO[[i]], names(smuGO)[i])
  return(eachMat)
}
GOMat <- as.data.frame(GOMat)
GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE)
GOTestWithCat <- GOTestWithCat[!is.na(GOTestWithCat$ontology), ]

## add ablog2FC
goSub <- smuGO[match(GOTestWithCat[, 1], names(smuGO))]
abLogFC <- sapply(goSub, function(x) {
  eachFC <- res[match(x, res$GeneID), 'log2FoldChange']
  return(mean(abs(eachFC), na.rm = TRUE))
})
GOTestWithCat$abLogFC <- abLogFC

## deal with NA and select BP MF and CC
termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write.csv(GOTestWithCat[GOTestWithCat$ontology == i, ],
            paste0('degseq4h_FC2_', i, '_withcat.csv'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## deal path
pathAnno <- getKEGGPathAnno('smu')
pathAnno[, 2] <- sapply(pathAnno[, 2], function(x){
  eachLen <- nchar(x)
  x <- substr(x, 1, eachLen - 29)
  return(x)
})

KEGGMat <- foreach(i = 1:length(smuKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(smuKEGG[[i]], names(smuKEGG)[i])
  return(eachMat)
}
KEGGMat <- as.data.frame(KEGGMat)

KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE)
KEGGTestWithCat$term <- pathAnno[match(KEGGTestWithCat[, 'category'], pathAnno[, 1]), 2]
KEGGTestWithCat$term %<>%
  strsplit(., split = ' - ', fixed  = TRUE) %>%
  sapply(., `[`, 1)
KEGGTestWithCat$ontology <- 'KEGG'

goSub <- smuKEGG[match(KEGGTestWithCat[, 1], names(smuKEGG))]
abLogFC <- sapply(goSub, function(x) {
  eachFC <- res[match(x, res$GeneID), 'log2FoldChange']
  return(mean(abs(eachFC), na.rm = TRUE))
})
KEGGTestWithCat$abLogFC <- abLogFC

write.csv(KEGGTestWithCat, file = 'degseq4h_FC2_KEGG_withcat.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

########################geneset plot###########################
setwd('/extDisk1/RESEARCH/smuSeqSongYing/kallisto_results')

library('ggplot2')
library('RColorBrewer')
library('latex2exp')
library('magrittr')

fname <- dir(pattern = 'withcat')
fbasename <- fname %>%
  strsplit(., split = '.', fixed = TRUE) %>%
  sapply(., `[`, 1)


for (i in seq_along(fname)) {
##for (i in 1:2) {
  pMat <- read.csv(fname[i], row.names = 1, stringsAsFactor = FALSE)

  ## pvalue < 0.05
  plotpMat <- pMat[, c('term', 'numDEInCat', 'over_represented_pvalue', 'abLogFC')]
  plotpMat[, 3] <- -log10(plotpMat[, 3])
  colnames(plotpMat) <- c('Name', 'Size', 'logpvalue', 'ablogFC')

  colorPal <- colorRampPalette(rev(c('red', 'yellow', 'cyan', 'blue')), bias=1)(10)

  ggplot(plotpMat[plotpMat[, 'logpvalue'] >= -log10(0.05), ], aes(x = ablogFC, y = Name)) +
    geom_point(aes(size = Size, colour = logpvalue)) +
    scale_size_continuous(name = 'Number of significant genes', range = c(3,8)) +
    scale_colour_gradientn(name = '-log10(P-value)', limits=c(0, max(plotpMat[, 3])), colours = colorPal) +
    ylab('') +
    xlab(TeX('Average |$\\log_{2}$FC|')) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 13, face = 'bold'))
  ggsave(paste0(fbasename[i], '.pdf'), width = 13)
}
###############################################################
