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
