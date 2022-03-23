awk -v cover=$2 '{if ($4>cover)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' $1
