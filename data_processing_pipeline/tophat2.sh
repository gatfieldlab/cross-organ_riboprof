#!/bin/bash

LIMIT=15
mkdir ./tophat


cat $1 | xargs -n 1 -P $LIMIT -I {} bash -c "echo Processing {} ... && "\
"mkdir ./tophat/{}/ && "\
"tophat2 --transcriptome-index=/local/databases/mouse/tophat/Mmusculus.GRCm38.75.dna.ensembl_data/known -p 3  -o tophat/{}/ Mmusculus.GRCm38.75.dna.ensembl {}_21-60.fastq && "\
"echo Finished processing {}"

echo "Done"

exit 0
