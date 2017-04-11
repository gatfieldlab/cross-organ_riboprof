#!/bin/bash
#violeta 09/04/14

# Get sample BASE name from first argument or die
if [ $# -eq 0 ] ; then
    echo 'Sample BASE name is missing!'
    exit 1
else
    BASE=$1
fi

# TR or RP
TY=${BASE:0:2}

# Bowtie mapping eff/sensitivity parameters
if [ $2 ] ; then
    BOWTIE_PARAMS=${*:2}
else
    BOWTIE_PARAMS="-L 15"
fi

# Declare location of filtered samples.
#Note that at this point the program has moved to mapped folder.
FILTERED="../filtered_data/"$BASE".fastq"

# 1. Map filtered samples against mouse rRNA, keep unaligned fastq (--un); keep output sam file (-S), suppress sam records for unaligned reads (--no-unal); keep log (2>).
# If no need to keep sam file, remove --no-unal -S sam/$BASE.vs_rRNA.sam and instead -S /dev/null
echo "Mapping $BASE against mouse rRNA"
bowtie2 $BOWTIE_PARAMS -x mouse_rRNA -q $FILTERED --un fq/$BASE.nonmrRNA.fastq -S sam/$BASE.vs_mrRNA.sam --no-unal 2> logs/$BASE.vs_mrRNA.bw2.log

#2. Map nonmrRNA.fastq against human rRNA, keep unaligned fastq; keep sam file (if needed) without unaligned read records; keep log.
echo "Mapping $BASE against human rRNA"
bowtie2 $BOWTIE_PARAMS -x human_rRNA -q fq/$BASE.nonmrRNA.fastq --un fq/$BASE.nonmhrRNA.fastq -S sam/$BASE.nonmrRNA_vs_hrRNA.sam --no-unal 2> logs/$BASE.non_mrRNA_vs_hrRNA.bw2.log

#3. Map nonmhRNA.fastq against tRNA, keep unaligned; keep sam without unaligned reads, keep log.
echo "Mapping $BASE against tRNA"
bowtie2 $BOWTIE_PARAMS -x Mmusculus_v38_tRNA -q fq/$BASE.nonmhrRNA.fastq --un fq/$BASE.nonmhrRNA_nontRNA.fastq -S sam/$BASE.nonmhrRNA_vs_tRNA.sam --no-unal 2> logs/$BASE.nonmhrRNA_vs_tRNA.bw2.log

#4. Map nonmhRNA_nontRNA.fastq against cDNA, keep unaligned; keep sam, without unaligned read records; keep log.
echo "Mapping $BASE against cDNA"
bowtie2 $BOWTIE_PARAMS -x Mmusculus.GRCm38.75.cdna.ensembl -q fq/$BASE.nonmhrRNA_nontRNA.fastq --un fq/$BASE.nonrRNA_noncDNA.fastq -S sam/$BASE.nonrRNA_vs_cDNA.sam --no-unal 2> logs/$BASE.nonrRNA_vs_cDNA.bw2.log

#5. Map nonrRNA_noncDNA.fastq against genome, keep unaligned; keep sam; keep log
echo "Mapping $BASE against genome"
bowtie2 $BOWTIE_PARAMS -x Mmusculus.GRCm38.75.dna.ensembl -q fq/$BASE.nonrRNA_noncDNA.fastq --un fq/$BASE.unmapped.fastq -S sam/$BASE.nonrRNA_noncDNA_vs_genome.sam 2> logs/$BASE.nonrRNA_noncDNA_vs_genome.bw2.log

#6. Remove fastq files except last unmapped.fastq
rm fq/$BASE.nonmrRNA.fastq
rm fq/$BASE.nonmhrRNA.fastq
rm fq/$BASE.nonmhrRNA_nontRNA.fastq
rm fq/$BASE.nonrRNA_noncDNA.fastq

echo "..finished"
exit 0
