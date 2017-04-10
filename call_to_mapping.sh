#!/bin/bash

#violeta 09/04/14

#Run in /filtered folder
# Run over filtered_files.txt containing fastq files names (of filtered samples)

#To do BEFORE:
	#Filter samples
	#Create /mapped folder

#This script calls bw2.mapping.sh. Both scripts saved in /local/violeta/scripts. If location changes, it has to be updated in this script.

#BOWTIE_PARAMS="-p 2 --very-sensitive --score-min L,0,-0.07 -k 20"
BOWTIE_PARAMS="-p 5 -L 15 -k 20"
LIMIT=9

#Copy the filtered_data.txt in /mapped and move to that folder. There, create sam,logs,fq folders.
cp $1 ../mapped
cd ../mapped
mkdir fq
mkdir logs
mkdir sam
cat $1 | xargs -n 1 -P $LIMIT -I {} bash -c "echo Processing {} ... && /local/violeta/scripts/bw2_mapping.sh {} $BOWTIE_PARAMS && "\
"echo Finished processing {}"

echo "Done."
exit 0

