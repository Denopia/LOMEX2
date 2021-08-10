#!/bin/bash

KMERS="data/sorted-oriented-readable-kmers-n3-s30-l5.txt"
READS="data/sl-40x-50k-i5-d5-s5.fastq"
OUTPUT="output/output.txt"
MN=3
ML=5
MW=30
K=$(($ML * $MN)) # ML*MN
A=5

/usr/bin/time -v ./lomex2  -k $KMERS -r $READS -o $OUTPUT -n $MN -l $ML -s $MW

