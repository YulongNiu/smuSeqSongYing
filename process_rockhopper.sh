SAMTOOLS_PATH=/home/Yulong/Biotools/samtools_1.8/bin/
BAM_PATH=/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/

cd ${BAM_PATH}

samfiles=($(ls | grep sam))

for j in "${samfiles[@]}"
do
    echo ${j}
    echo ${j%%.*}.bam
    ${SAMTOOLS_PATH}/samtools sort -@ 8 -o ${j%%.*}.bam ${j}
    ${SAMTOOLS_PATH}/samtools index ${j%%.*}.bam
    rm ${j}
done
