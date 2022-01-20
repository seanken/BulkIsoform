transfa=$1 
fa=$2
name=$3


echo Add on decoys
cat $transfa $fa > comb.fa
grep "^>" $fa | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt

echo Make Reference
salmon index -t comb.fa -i $name --decoys decoys.txt
rm comb.fa
