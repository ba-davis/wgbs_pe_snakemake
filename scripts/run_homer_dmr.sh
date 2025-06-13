#!/bin/bash

inpath=$1
genome=$2
ens_ref=$3
outdir_base=$4

# if ens_ref is true, remove the chr from the input bed files
# sed command

for file in $inpath/*.GREAT.bed; do
    echo $file
    name=${file##*/}
    echo $name
    outdir=$(basename "$file" .GREAT.bed)
    echo $outdir
    outdir_full="${outdir_base}${outdir}"
    echo $outdir_full
    mkdir $outdir_full

    # if enembl ref genome, remove chr from input bed file
    if [ "$ens_ref" = "TRUE" ]; then
	echo "Ensembl ref genome is set as TRUE"
	echo "Removing chr string from input bed file"
	sed 's/^chr\([0-9XYM]*\)/\1/' $file > ${file/.GREAT./.HOMER.}
	file=${file/.GREAT./.HOMER.}
    else
	echo "Reference genome is not Ensembl"
    fi

    # run homer
    findMotifsGenome.pl $file $genome $outdir_full -size given -len 8,10,12 -p 8
done

echo "homer complete." > ${outdir_base}homer_complete.txt
