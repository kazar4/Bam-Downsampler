# Bam-Downsampler
Script to preform downsampling on a BAM file

## Usage

`python3 downsample_bam.py`

### Required attributes
-i filePathToBamFiles (can use * notation to downsample multiple at once, can give sam or bam files)

-d Amount (the amount of each primer to be left after downsampling)

-p filePathToPrimers (text file containing the primers to keep, example primers file at end)

### Optional attributes
-s (sorts the bam file after downsampling)

## Example:
`python3 downsample_bam.py -i HLA-LA2/*.bam -d 10 -p primers.txt -s`

## Primers file format example
#### Can have up to two forward and reverse strands, seperated with a comma

```
TARGET	FORWARD	REVERSE
HLA-A	ATCCTGGATACTCACGACGCGGAC	CATCAACCTCTCATGGCAAGAATTT
HLA-B	AGGTGAATGGCTCTGAAAATTTGTCTC,TAG	AGAGTTTAATTGTAATGCTGTTTTGACACA,GAT
HLA-C	CAGCACGAAGATCACTGGAA	TGAGGAAAAGGAGCAGAGGA
```

## Requirements
Command line installation of `samtools`

pandas, argparse, csv

Script should work on any system with python as long as samtools is installed, ubuntu directions for samtools below

`apt-get install samtools`
