# uniqueKmer
identification of kmers with n occurrences from metagenomic databases

## Installation
make

## Sample command to run
./uniqueLoci -i <directory containing > -o <output file name> -k <kmer size> -e 0.0000001 -d 0.0005 -t <number of threads> -c <maximum kmer count to report in the output>

## Warnings

1. Each thread will create it's own copy of the sketch, meaning the total memory occupied by the program is size of single sketch (see below) times number of threads

Sample epsilon and delta parameters and total size of a sketch
### Note: sensitivity will depend on the size of the input data, so setting these parameters requires some intuition about the sketch.
### Lower values of epsilon and delta will result in higher sensitivity at higher memory (epsilon) or compute time (delta)
| Epsilon   | Delta | Total Sketch Size (MB) |
| --------- |:-----:|:----------------------:|
| 0.0001 (27183 elements)   | 0.005(6 hashes) | 0.3 |
| 0.00001 (271829 elements)  | 0.005(6 hashes) | 3 |
| 0.000001 (2718282 elements) | 0.005(6 hashes) | 31 |
| 0.0000001 (27182818 elements)| 0.005(6 hashes) | 311 |
