#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: test_sample_seq_prep.sh

declare -A plates_per_country
declare -A controls_per_plate

first_line=true

while IFS=, read -r Sample_ID Plate Origin

do
    if [ "$first_line" = true ] ; then
        first_line=false
        # echo "Skipping first line..."
        continue
    fi 

    if [[ "$Origin" == "NTC" || "$Origin" == "Mock" || "$Origin" == "PosC" ]] ; then
        # echo "Origin is $Origin, skipping..."
        controls_per_plate[$Plate]="${controls_per_plate[$Plate]}
        $Sample_ID"
        continue

    fi  
    # echo "Debug: $Origin ; $Sample_ID ; $Plate"
    mkdir -p $Origin/
    mv $Sample_ID/ $Origin/
    plates_per_country[$Origin]="${plates_per_country[$Origin]}
    $Plate"
    # delete if duplicate
    plates_per_country[$Origin]="$(echo "${plates_per_country[$Origin]}" | sort -u)"

done < <(csvcut -c Sample_ID,Plate,Origin $1.csv)

echo "${!plates_per_country[@]}"
echo "${plates_per_country[@]}"

for country in ${!plates_per_country[@]}; do
	plates="${plates_per_country[$country]}"
	echo "$country -> ${plates}"
    for plate in $plates; do
        controls="${controls_per_plate[$plate]}"
        cp -r $controls $country/
    done
done

