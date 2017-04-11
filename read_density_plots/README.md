# DENSITY PLOTS #
Scripts to generate transcript ribosome and RNA-seq read density plots.

### process_reads_v7.py ###
Run with option -c 5prime this script counts reads per position (respect to transcript's 5'end)
Output needs to be sorted by gene (using qsort) and indexed with: tr_index_py3_mod.py

### get_density4.py ###
This script extracts the density information of a selected transcript id (input is a sorted 5-prime count file).
This script is necessary for plot_density_v4_kidney.py

### plot_density_v4_kidney.py ###
Script to generate read density plots. Takes one transcript ID.
With options to draw manual boxes (fir example, for uORFs)
