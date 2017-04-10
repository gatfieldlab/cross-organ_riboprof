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

circadian analysis
babel
normalization, rpkm and te?
hellinger distance
euclidean distance and clustering

general functions.R

