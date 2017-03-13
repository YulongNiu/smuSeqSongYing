#!/bin/sh

INDEX_BASE=/home/Yulong/Biotools/RefData/smu/
NAME_PREFIX=NC_004350
INDEX_NAME=${NAME_PREFIX}_kallisto_index

cd ${INDEX_BASE}

## build index
kallisto index --index ${INDEX_NAME} \
         ${NAME_PREFIX}_transcripts.fa 
