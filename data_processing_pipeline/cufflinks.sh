#!/bin/bash

LIMIT=9
#mkdir ./tophat
mkdir ./cufflinks


cat $1 | xargs -n 1 -P $LIMIT -I {} bash -c "echo Processing {} ... && "\
"mkdir ./cufflinks/{}/ && "\
"cufflinks -o cufflinks/{}/ --GTF /local/databases/mouse/gtf/Mus_musculus.GRCm38.75.gtf --frag-len-mean 37 --frag-len-std-dev 8 --compatible-hits-norm --multi-read-correct --upper-quartile-norm --frag-bias-correct /local/databases/mouse/fasta/Mmusculus.GRCm38.75.dna.ensembl.fa -p 3 tophat/{}/accepted_hits.bam && "\
"echo Finished processing {}"

echo "Done"

exit 0
