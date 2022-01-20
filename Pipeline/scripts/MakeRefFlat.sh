gtf=$1
output=$2

gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $gtf temp.txt
awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'  temp.txt > $output
