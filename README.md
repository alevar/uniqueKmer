# uniqueKmer
identification of kmers with n occurrences from metagenomic databases

## Installation
make

## Sample command to run
`./uniqueLoci -i <directory containing metagenomic samples> -o <output file name > -k <kmer size> -e 0.00000001 -d 0.0005 -t <number of threads> -c <maximum kmer count to report in the output>`

## Warnings

1. Each thread will create it's own copy of the sketch, meaning the total memory occupied by the program is: size of single sketch (see below) times number of threads
2. No ambiguous characters are currently supported: only ACTGN alphabet. Kmer's containing ambiguous bases will not be counted.
3. reverse complement is automatically considered for each kmer and unique occurrence means that there is no other kmer with a reverse complement identical to current kmer
4. Input directory must have the following structure: each reference 
5. Each thread will output its own comma-separated file with a threadID appended to the basename (output argument). User may wish to concatenate the files together into a single file using `cat file1 file2 ... fileN > result`. Full description of the file format is provided below

## Output format - CSV
1. file path to the genome
2. fasta sequence name within the file
3. 0-based coordinate of the first nucleotide in the kmer
4. number of occurrences of the kmer
5. kmer sequence

## Sample epsilon and delta parameters and total size of a sketch
### Note: sensitivity will depend on the size of the input data, so setting these parameters requires some intuition about the sketch.
### Lower values of epsilon and delta will result in higher sensitivity at higher memory (epsilon) or compute time (delta)
| Epsilon   | Delta | Total Sketch Size (MB) |
| --------- |:-----:|:----------------------:|
| 0.0001 (27183 elements)   | 0.005(6 hashes) | 0.2 |
| 0.00001 (271829 elements)  | 0.005(6 hashes) | 1.5 |
| 0.000001 (2718282 elements) | 0.005(6 hashes) | 16 |
| 0.0000001 (27182818 elements)| 0.005(6 hashes) | 155 |
