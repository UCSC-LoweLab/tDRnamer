#!/usr/bin/env bash

REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )

# Download example data
echo -e "\nDownloading example data\n"
wget "http://trna.ucsc.edu/tDRnamer/data/examples/GM12878_AlkB_rep1_merged.fastq.gz"
wget "http://trna.ucsc.edu/tDRnamer/data/examples/ExampleSequences.fa"
wget "http://trna.ucsc.edu/tDRnamer/data/examples/ExampleNames.txt"

# Download hg38 reference genome sequence
echo -e "\nDownloading reference genome\n"
wget "http://trna.ucsc.edu/tDRnamer/data/examples/hg38_primary_assembly.fa.gz"
gzip -d hg38_primary_assembly.fa.gz

# Get tRNA information
echo -e "\nDownloading tRNA data\n"
wget "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"
tar xvf hg38-tRNAs.tar.gz

# Create reference database
echo -e "\nCreating reference database\n"
mkdir -p db
"$SCRIPTDIR/create_tDRnamer_db" -d db/hg38 -g hg38_primary_assembly.fa -t hg38-tRNAs-detailed.out -s hg38-tRNAs-detailed.ss -n hg38-tRNAs_name_map.txt -q

# Name tDRs from fastq
echo -e "\nExample Run 1\n"
mkdir -p run1
"$SCRIPTDIR/tDRnamer" -m seq -s GM12878_AlkB_rep1_merged.fastq.gz -d db/hg38 -r euk -o run1/GM12878_AlkB_rep1

# Name tDRs from fasta
echo -e "\nExample Run 2\n"
mkdir -p run2
"$SCRIPTDIR/tDRnamer" -m seq -s ExampleSequences.fa -d db/hg38 -r euk -o run2/ExampleSequences

# Name tDRs from fasta with maximum sensitivity
echo -e "\nExample Run 3\n"
mkdir -p run3
"$SCRIPTDIR/tDRnamer" -m seq -s ExampleSequences.fa -d db/hg38 -r euk --max -o run3/ExampleSequences

# Name tDRs from fasta with nucleotide variations                                                                                                                                                                    
echo -e "\nExample Run 4\n"
mkdir -p run4
"$SCRIPTDIR/tDRnamer" -m seq -s ExampleSequences.fa -d db/hg38 -r euk --var -o run4/ExampleSequences

# Find tDRs by names
echo -e "\nExample Run 5\n"
mkdir -p run5
"$SCRIPTDIR/tDRnamer" -m name -n ExampleNames.txt -d db/hg38 -r euk -o run5/ExampleNames
