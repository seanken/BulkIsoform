fa=$1
gtf=$2
name=$3

rsem-prepare-reference --gtf $gtf --star $fa $name
