# DATA ANALYSIS AND PLOTTING #

### functions.R
Script containing functions used throughout R analyses.

### normalisation_rpkm_te.R
Read count files, filtering, upper-quartile normalisation, pseudo-rpkm and Translation efficiency calculations.

## CIRCADIAN ANALYSES AND BABEL ##
### circadian_analysis_kidney.R
Script used for rhythmicity detection in RPF and RNA datasets, daily profile plot generation transcriptome-wide, and Babel analysis for translation efficiency regulation.

## TRI-NUCLEOTIDE PERIODICITY AND FRAME BIAS ##
### calc_codons_fractions.py
This script calculates the fraction of reads mapping to each nucleotide of the specified region, averaged over all input files. Input files are 5-prime sorted count files. The regions used in the study were [-20:+200] around start codon (specified as -r "*cds-20:*cds+200") and [-200:+20] around stop codon (-r "cds*-200:cds*+20") for tri-nucleotide periodicity analyses, and whole cds (-r "*cds:cds*") for frame bias analysis.

## uORF PREDICTION ##
### uorf_multiprocess_kidney_v2.py
This script takes a list of transcript ids and predicts the presence of AUG-starting uORFs in 5'UTRs. Outputs the start and stop of potential uORFs, the coverage, and the normalized read counts in each frame.

### translated_uorfs.R
This script predicts translated uORFs, based on a  minimun coverage of 10 %, frame bias, and preferred first frame (respect to uORF start).

## TRANSCRIPT USAGE DIVERSITY ##
Hellinger distance was used to quantify the dissimilarity in relative transcript isoform expression between organs, which was later correlated with differences in translation efficiency across tissues.
### get_fmi.py
This script gathers the relative expression proportions of all expressed protein-coding transcripts (estimated from our RNA-seq, see cufflinks.sh).
### hellinger.R
This script reads the output of get_fmi.py and calculated hellinger distances across organs for each gene.
### transcript_diversity_v3.py
In order to detect the transcript features that were associated with tissue specificity in TE, we selected genes whose transcript diversity between both organs originated from, or was excluded from, 5′ UTR, CDS, or 3′ UTR, based on feature annotation information for the detected protein-coding transcripts.
### transcript_features_v2.py
For single-isoform genes, we investigated whether a particular transcript characteristic (length, GC content, Kozak context, structure) could be predictive of differential TE.

## EUCLIDEAN DISTANCES AND HIERARCHICAL CLUSTERING  OF RHYTHMIC GENES ##
### euclidean_clustering.R
Script with functions to calculate dissimilarity matrices between the four datasets (RNA and RPF, liver and kidney) for genes of interest, draw clustering trees for individual genes, and a weighted average dissimilarity matrix (and tree) for several genes.
