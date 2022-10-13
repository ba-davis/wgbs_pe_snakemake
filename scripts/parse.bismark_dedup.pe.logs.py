#---------------------------------------------------------------------------------#
# Purpose: to parse relevant info from a directory of bismark dedup log files     #
#          and produce an outfile containing a table of stats per sample          #
#                                                                                 #
# Useage: ./parse.bismark_dedup.pe.logs.py -d /path/to/reports -o path/to/outfile #
#---------------------------------------------------------------------------------#

# Notes:
#       - Assumes the bismark log files end with "report.txt"
#       - Written for PE files (looks for "Sequence pairs analysed in total" line)
#       - Don't need path to trimmed files anymore

# TO DO: - Make it PE or SE agnostic? Grab for something that is in both reports
#       - OR make a PE version and a SE version


# Import libraries
import optparse
import pandas as pd
import re
import os


# Initialize command line options for input directory and output file
p=optparse.OptionParser()
p.add_option("-d", action = "store", dest = "directory")
p.add_option("-o", action = "store", dest = "outfile")

opts,args=p.parse_args()
in_dir=opts.directory
outfile=opts.outfile

# open outfile for writing
fhw = open(outfile, "w+")
# write colnames to outfile
fhw.write('Sample_Name' + '\t' + 'Input_Alignments' + '\t' + 'Dup_Aln_Rm' + '\t' + 'Dup_Aln_Rm_Perc' + '\t' + 'Dup_Aln_Positions' + '\t' + 'Deduped_Alns_Remaining' + '\t' + 'Deduped_Alns_Remaining_Perc' + '\n')

# Loop through files in input directory, and work on files ending in *report.txt
#   to grab lines of interest and record the metric to an output file
for file in sorted(os.listdir(in_dir)):
    filename = os.fsdecode(file)
    if filename.endswith("report.txt"):
        path_file=in_dir + '/' + filename
        theFile = open(path_file,'r')
        FILE = theFile.readlines()
        theFile.close()
        printList = []
        for line in FILE:
            # Obtain sample name from bam file)
            if ('alignments analysed in' in line):
                intname=line.split()[6]
                intname2=intname.split('/')[-1]
                intname3=re.sub("(_val_1.*)$", "", intname2)
                printList.append(intname3 + '\t')
            # Obtain Number of Alignments Analyzed
            if ('alignments analysed in' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Obtain Number of Dup Alignments removed
            if ('duplicated alignments removed' in line):
                intname=line.split('\t')[1].rstrip()
                intname2=intname.split(' ')[0]
                printList.append(intname2 + '\t')
            # Obtain Percentage of Dup Alignments removed
            if ('duplicated alignments removed' in line):
                intname=line.split('\t')[1].rstrip()
                intname2=intname.split(' ')[1]
                intname3=intname2.replace('(','')
                intname4=intname3.replace(')', '')
                printList.append(intname4 + '\t')
            # Obtain Number of positions dup alignments occurred at
            if ('were found at' in line):
                intname=line.split('\t')[1]
                intname2=intname.split(' ')[0]
                printList.append(intname2 + '\t')
            # Obtain Number of dedup alignments remaining
            if ('leftover sequences' in line):
                intname=line.split(': ')[1]
                intname2=intname.split(' ')[0]
                printList.append(intname2 + '\t')
            # Obtain Percent of dedup alignments remaining
            if ('leftover sequences' in line):
                intname=line.split(': ')[1]
                intname2=intname.split(' ')[1]
                intname3=intname2.replace('(','')
                intname4=intname3.replace(')', '')
                printList.append(intname4 + '\n')
        for item in printList:
            fhw.write(item)

    else:
        continue

fhw.close()
      