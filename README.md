# routines-with-sequences

streamlining_aln.sh 
# (1) What it does and why
# This script performs multiple sequence alignment using muscle with default parameters.
# It assumes that you have 'muscle' executable in your path.
# Running this script is nothing special except that it will first dereplicates the input, align the unique sequences, then re-replicates.
# I found that using it reduces running time substantially, compared to running muscle on the whole input sequences.
# But that is surely because my data sets (i.e. intra-specific bacterial isolate collection) often contains identical sequences.
# There will be no benefit if your input data is already non-redundant (i.e. one sequence per species, or per group defined at certain level)
#
# (2) DereplicateRereplicateSeqFasta.py
# The shell script (streamlining_aln.sh) uses this python script (DereplicateRereplicateSeqFasta.py)
# This rudimentary script was made by me who began to learn bioinformatics coding.
#
# (3) Usage
./streamlining_aln.sh [Input fasta path] [Output fasta path]
#
# (4) Warning
# There is no help message. Don't try -h --help or fancy things like that.
