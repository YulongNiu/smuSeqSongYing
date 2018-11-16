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

fmember=('4h_sm_1' '4h_sm_2' '4h_sm_3' '4h_srta_1' '4h_srta_2' '4h_srta_3' '24h_sm_1' '24h_sm_2' '24h_sm_3' '24h_srta_1' '24h_srta_2' '24h_srta_3')

${KALLISTO_PATH}/kallisto quant \
                -t 8 \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                -o kallisto_results/4h_sm_1 \
                cleandata/4h_sm_1_R1.fq.gz clearndata/4h_sm_1_R2.fq.gz

${KALLISTO_PATH}/kallisto quant \
                -t 8 \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                -o kallisto_results/24h_sm_1 \
                cleandata/4h_sm_1_R1.fq.gz cleandata/24h_sm_1_R2.fq.gz

date
