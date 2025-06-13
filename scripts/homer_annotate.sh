#!/bin/bash

inpath=$1
genome=$2
ens_ref=$3
outdir_base=$4

# first, concatenate all the .motif fules in the knownResults directory into one .motif file
# cat *.motif > knownResults.all.motif

# Use a for loop to iterate over all 'knownResults' directories
for dir in ${outdir_base}*/knownResults; do
    # Concatenate all .motif files in the current knownResults directory into knownResults.all.motif
    echo $dir
    cat "$dir"/*.motif > "$dir"/knownResults.all.motif
done

# loop through input files
for file in $inpath/*.GREAT.bed; do
    outdir=$(basename "$file" .GREAT.bed)
    outdir_full="${outdir_base}${outdir}"
    #echo $outdir
    #echo $outdir_full

    # if enembl ref genome, look for the HOMER files instead of GREAT
    if [ "$ens_ref" = "TRUE" ]; then
	echo "Ensembl ref genome is set as TRUE"
	echo "Looking for HOMER input file"
	file=${file/.GREAT./.HOMER.}
    else
	echo "Reference genome is not Ensembl"
    fi

    # define motif file
    motif_file=$outdir_full/knownResults/knownResults.all.motif

    # run homer annotate
    annotatePeaks.pl $file $genome -size given -m $motif_file > ${outdir_base}$outdir.knownMotifs.all.annot.output.txt
    
done

echo "homer annotate complete." > ${outdir_base}homer_annotate_complete.txt
