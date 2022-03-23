awk -v cover=$2 '{if ($3 - $2 > cover) print }' $1 | awk '{print $s "\tnovel"NR"_"($3-$2+1)"bp"}'

