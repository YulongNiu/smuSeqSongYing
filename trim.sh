date

BASE_PATH=/extDisk1/RESEARCH/smuSeqSongYing
RAW_PATH=${BASE_PATH}/rawdata
CLEAN_PATH=${BASE_PATH}/cleandata

FASTP_PATH=/home/Yulong/Biotools/fastp

CORENUM=8

fmember=('4h_sm_1' '4h_sm_2' '4h_sm_3' '4h_srta_1' '4h_srta_2' '4h_srta_3' '24h_sm_1' '24h_sm_2' '24h_sm_3' '24h_srta_1' '24h_srta_2' '24h_srta_3')

cd ${RAW_PATH}

for i in "${fmember[@]}"
do
    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -g -p -c \
                 -h ${i}.html -j ${i}.json \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz

done

date
