# DATA PROCESSING PIPELINE #

## ADAPTOR TRIMMING ## 
### remove_adaptors2.sh
Sequencing adaptors were removed with cutadapt v1.3 with the following options:
-a AGATCGGAAGAGCACACGTCT --match-read-wildcards -m 6

## SIZE FILTERING ## 
### size_filtering.py
RPF samples were filtered to keep reads of lengths [26:35], and RNA samples for read size lengths [21:60].

## MAPPING ##
Two mapping strategies were applied: a sequential mapping to obtain transcriptome alignments; and a direct genome mapping to estimate expressed transcript models from our RNA-seq data.

### SEQUENTIAL MAPPING ### 
### bw2_mapping.sh
Trimmed and filtered reads were sequentially mapped to: 1) Mouse rRNA, 2) Human rRNA, 3) Mouse tRNA,4) Mouse cDNA (Ensembl release 75) and 5) Mouse genome  (Genome Reference Consortium GRCm38.p2).
After each alignment, only non-mapping reads were used for the following mapping. Only mouse cDNA-mapping sequences were used for subsequent analysis.
### filter_sam.v2.py
When multiple alignment occured, only alignments with maximun alignment scores (AS) were kept.

### GENOMIC MAPPING ##
### tophat2.sh
Trimmed and filtered reads were mapped directly agains the mouse genome with Tophat2.
### cufflinks.sh
Cufflinks was then used to estimate the number of fragments per kilo base of exon per million fragments map (FPKM) for each transcript (of Ensembl release 75).
### parse_cuffcmp.py
Resulting expression estimates were filtered to keep transcript isoforms with an FPKM > 0.1, a LOW 95 > 0.05, and a FMI > 2.0 in at least 3 samples.
### prepare_gtf_v3.py
A database of expressed transcripts/genes was built based on these estimations and used in all subsequent read counting and analyses.

## COUNTING ##
### process_reads_v6.py
Reads were counted towards their location within annotation feature (5'UTR, CDS, 3'UTR) per gene. 
Only reads that mapped uniquely to a single expressed gene were counted. 
For reads mapping to multiple transcript isoforms of the same gene (multiple-isoform genes), reads that did not map unambiguosly to a single feature were assigned in the following preference order: CDS > 5'UTR > 3'UTR.

This script was also used to count reads with positional information with respect to the transcript's 5' end (with the options -c 5prime).
