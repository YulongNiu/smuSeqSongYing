date

BEDTOOLS_PATH=/home/Yulong/Biotools/bedtools2/bin

INDEX_PATH=/home/Yulong/Biotools/RefData/smu

cd ${INDEX_PATH}

${BEDTOOLS_PATH}/bedtools getfasta -fi NC_004350.fna -bed NC_004350.gff -name -s -fullHeader -fo NC_004350_cdna.fa

date
