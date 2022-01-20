gtf=$1
outfile=$2


grep -v \# $gtf | awk '{if($3=="transcript"){print $0}}'  | awk -F 'gene_name ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g'> gene.txt
grep -v \# $gtf | awk '{if($3=="transcript"){print $0}}'  | awk -F 'transcript_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g'> trans.txt
paste trans.txt gene.txt | sort | uniq > $outfile
rm gene.txt
rm trans.txt
