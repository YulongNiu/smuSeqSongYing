##################################proprocess annotation###############
##~~~~~~~~~~~~~~~~~~~~~~~~~~~gff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smugff <- read.delim('/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/genomes/Streptococcus_mutans_UA159/NC_004350.gff', skip = 3, header = FALSE, stringsAsFactor = FALSE)

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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~add ptt/rnt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



######################################################################
