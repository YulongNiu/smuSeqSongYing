date

KALLISTO_PATH=/home/Yulong/Biotools/kallisto_linux-v0.44.0

INDEX_PATH=/home/Yulong/Biotools/RefData/smu

PROJECT_PATH=/extDisk1/RESEARCH/smuSeqSongYing

# ## build index
# ${KALLISTO_PATH}/kallisto index \
#                 -i ${INDEX_PATH}/smu_kallisto_idx \
#                 ${INDEX_PATH}/Streptococcus_mutans_ua159.ASM746v2.cdna.all.fa.gz


## quantification
cd ${PROJECT_PATH}

${KALLISTO_PATH}/kallisto quant \
                -t 8 \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                -o kallisto_results/4h_sm_1 \
                rawdata/4h_sm_1_R1.fq.gz rawdata/4h_sm_1_R2.fq.gz

${KALLISTO_PATH}/kallisto quant \
                -t 8 \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                -o kallisto_results/24h_sm_1 \
                rawdata/24h_sm_1_R1.fq.gz rawdata/24h_sm_1_R2.fq.gz

date
