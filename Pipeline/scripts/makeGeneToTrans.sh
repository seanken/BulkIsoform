gtf=$1
gene_name=$2
outfile=$3

if [[ "$gene_name" == "gene_name" ]]
then
grep -v \# $gtf | awk '{if($3=="transcript"){print $0}}'  | awk -F 'gene_name ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g'> gene.txt
else
grep -v \# $gtf | awk '{if($3=="transcript"){print $0}}'  | awk -F 'gene_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g'> gene.txt
fi
grep -v \# $gtf | awk '{if($3=="transcript"){print $0}}'  | awk -F 'transcript_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g'> trans.txt
paste trans.txt gene.txt | sort | uniq > $outfile
rm gene.txt
rm trans.txt
