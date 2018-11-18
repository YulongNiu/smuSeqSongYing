date

KALLISTO_PATH=/home/Yulong/Biotools/kallisto_linux-v0.44.0

INDEX_PATH=/home/Yulong/Biotools/RefData/smu

PROJECT_PATH=/extDisk1/RESEARCH/smuSeqSongYing
CLEAN_PATH=${PROJECT_PATH}/cleandata
RES_PATH=${PROJECT_PATH}/kallisto_results

## build index
${KALLISTO_PATH}/kallisto index \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                ${INDEX_PATH}/NC_004350_cdna_name.fa

## quantification
cd ${CLEAN_PATH}


fmember=('4h_sm_1' '4h_sm_2' '4h_sm_3' '4h_srta_1' '4h_srta_2' '4h_srta_3' '24h_sm_1' '24h_sm_2' '24h_sm_3' '24h_srta_1' '24h_srta_2' '24h_srta_3')

for i in "${fmember[@]}"
do
    echo "Quantify ${i}"
    ${KALLISTO_PATH}/kallisto quant \
                -t 8 \
                -i ${INDEX_PATH}/smu_kallisto_idx \
                -o ${RES_PATH}/${i} \
                ${i}_R1.fq.gz ${i}_R2.fq.gz
done

date
