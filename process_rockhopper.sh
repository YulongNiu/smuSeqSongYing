SAMTOOLS_PATH=/home/Yulong/Biotools/samtools_1.8/bin/
BAM_PATH=/extDisk1/RESEARCH/smuSeqSongYing/Rockhopper_Results/

cd ${BAM_PATH}

bamfiles=($(ls | grep sam))

for j in "${bamfiles[@]}"
do
    echo ${j}
    echo ${j%%.*}.bam
    ${SAMTOOLS_PATH}/samtools sort -@ 8 -o ${j%.*}.bam ${j}
done
