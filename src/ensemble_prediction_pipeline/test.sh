#!/bin/bash
# INPUT='data/processed/ensemble-sparkle-params-regress-rnai-matched.csv'
# OLDIFS=$IFS
# IFS=','
# [ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
# {
#   read
#   while read start end model
#   do
#   	echo "Start : $start"
#   	echo "End : $end"
#   	echo "Model : $model"
#   done < $INPUT
#   IFS=$OLDIFS
# } < $INPUT

#!/bin/bash
filename="data/processed/ensemble-sparkle-params-regress-rnai-matched.csv"
{
    read
    while IFS=, read -r start end model
    do
      echo "Start : $start"
      echo "End : $end"
      echo "Model : $model"
    done
} < $filename
