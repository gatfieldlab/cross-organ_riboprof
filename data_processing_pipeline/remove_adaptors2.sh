#!/bin/bash


LIMIT=15
TRIMMED="../trimmed_data"

cat $1 | xargs -n 1 -P $LIMIT -I {} bash -c "echo Processing {} ... && gzip -d -c {}.fastq.gz | "\
"cutadapt -a AGATCGGAAGAGCACACGTCT --match-read-wildcards -m 6 "\
"-o $TRIMMED/{}_trimmed.fastq -y :trimmed "\
"- > $TRIMMED/{}_trimming_summary.txt && "\
"echo Finished processing {}"

echo " Script completed. "
exit 0
