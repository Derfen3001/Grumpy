# Grumpy
HiC data processing

This script takes as imput HiC data in .cool format and a file containing the genomic regins of interest (for example TADs).
As result the program calculates the interactions between and within the provided genomic regions.

#### List of genomic regions

This is a tsv without header containing chr, start, end, Type_1 and Type_2.
Type_1 and Type_2 are mandatory colums that identify characteristics of the genomic regions (for exampem Type_1: Hihg Low H3K9me2 and Type_2: High Low H3k27me3).
If not needed the insert a unique cluster each column.

#### Usage

grumpy.py [-h] -i INPUT [INPUT ...] -l LIST [-p PROCESSORS]
                 [-o OUTPUTDIR] [-f ELEMENT] [-b] [-pl] [-s] [-X] [-coo]
                 [-trans] [--version]
             
  -h, --help            show this help message and exit
  
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]\n
                        One or multiple cooler files
  
  -l LIST, --list LIST  Sorted TAB list with no header and Chr Start End
                        Type_1_clusters Type_2_clusters colums
  
  -p PROCESSORS, --processors PROCESSORS
                        number of cores to use, if not specified it uses all
                        the processors available
  
  -o OUTPUTDIR, --output OUTPUTDIR
                        output directory
  
  -f ELEMENT, --feature_name ELEMENT
                        Name of the genomic feature
  
  -b, --boundaries      calculates interaction numbers over consecutive
                        elements (the List provided must be sorted)
  
  -pl, --plot           activates exploratory plotting options
  
  -s, --sex_Chr         Keep sexual chromosomes
  
  -X, --ChrX            Keep X chromosome
  
  -coo, --Coordinates   Recover Feature number with genomic positions
  
  -trans, --incl_trans  Include in trans interactions, False by default
  
  --version             show program's version number and exit
