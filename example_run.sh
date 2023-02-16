#!/bin/bash

# Example SAKE pipeline
# Tested on a Linux system

# Numnber of minimizer in a strobemer
N=3
# Minimizer length in a strobemer
L=7
# Minimizer window length in a strobemer
W=11
# Distance between the starting position of two consecutive minimizer windows
V=25
# Proportional bundling threshold
B=0.8
# Proportional character support threshold
Z=0.4
# Minimizer-k-mer set abundance threshold
A=4
# Bundle size threshold
G=4
# Character support threshold
C=4
# Number of threads for BFC
BT=1
# Number of threads for SAKE
LT=1
# Approximate number of k-mers, required by BFCounter
P=10000000
# Abundance for k-mer orienter, filters out k-mers with lesser abundance 
KA=1
# Strobemer length
SL=$((L * N))
# Default k-mer length, according to the formula in the manuscript (maximum size). Can be made smaller. If you make it larger, might cause problems.
K=$(($((N - 1))* V - W + L + 1))
# Maximum number of sequences allowed in POA building, Strobemers with greater number of occurrences are ignored
POA_LIMIT=10000000

# Path to reads
READS_PATH="example_data/reads40x_0.05error.fastq"
# Path to BFCounter binary strobemers
STROBEMERS_BIN_PATH="example_output/bfcounter_strobemers.bin"
# Path to BFCounter readable strobemers
STROBEMERS_PATH="example_output/bfcounter_strobemers.txt"
SAKE_OUTPUT_PATH="example_output/sake_sequences.txt"
SAKE_OUTPUT_AS_KMERS="example_output/sake_"$K"-mers.txt"
SAKE_OUTPUT_AS_KMERS_ORIENTED="example_output/sake_"$K"-mers_oriented.txt"
SAKE_OUTPUT_AS_KMERS_ORIENTED_SORTED="example_output/sake_"$K"-mers_oriented_sorted.txt"
SAKE_OUTPUT_AS_KMERS_ORIENTED_SORTED_UNIQ="example_output/FINAL-"$K"-mers.txt"

echo "k-mer length is"
echo $K

# Run BFCounter count
/usr/bin/time -v ./BFCounterForStrobemers/BFCounter count --verbose -k $SL -n $P -t $BT -m $N -l $L -w $W -f $V -o $STROBEMERS_BIN_PATH $READS_PATH
# Run BFCounter dump
/usr/bin/time -v ./BFCounterForStrobemers/BFCounter dump -k $SL -i $STROBEMERS_BIN_PATH -o $STROBEMERS_PATH
# Run SAKE
/usr/bin/time -v ./SAKEandSPOA/sake/sake -k $STROBEMERS_PATH -r $READS_PATH -o $SAKE_OUTPUT_PATH -n $N -l $L -s $W -f $V -b $B -z $Z -a $A -g $G -c $C -t $LT -p $POA_LIMIT
# Split SAKE output sequences into k-mers
/usr/bin/time -v python ./pytools/split_into_kmers.py $SAKE_OUTPUT_PATH $SAKE_OUTPUT_AS_KMERS $K
# Orient SAKE k-mers
/usr/bin/time -v python ./pytools/kmer_orienter.py $SAKE_OUTPUT_AS_KMERS $SAKE_OUTPUT_AS_KMERS_ORIENTED $K $KA
# Sort SAKE k-mers
/usr/bin/time -v sort $SAKE_OUTPUT_AS_KMERS_ORIENTED > $SAKE_OUTPUT_AS_KMERS_ORIENTED_SORTED
# Remove any duplicate k-mers
/usr/bin/time -v uniq $SAKE_OUTPUT_AS_KMERS_ORIENTED_SORTED > $SAKE_OUTPUT_AS_KMERS_ORIENTED_SORTED_UNIQ


