GTF_PATH=/home/Yulong/Biotools/RefData/smu/NC_004350.gff
MAP_PATH=/home/Yulong/RESEARCH/SongYing_MJ201409021010/Rockhopper_Results/

cd ${MAP_PATH}

featureCounts -T 8 -p -a ${GTF_PATH}\
              -t "Coding gene" \
              -g gene_id \
              -o counts.txt sm_1.sam



