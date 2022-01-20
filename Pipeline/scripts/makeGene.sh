gtf=$1
bed=$2

grep -v \# $gtf | awk -v OFS='\t' '{if($3=="gene"){print $1"\t"$4"\t"$5"\t"$14"\t.\t"$7}}' | sed 's/;//g' | sed 's/\"//g' > $bed

