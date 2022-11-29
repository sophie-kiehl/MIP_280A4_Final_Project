# MIP 280A4 Final Project 2022

This report documents the final project I did for MIP 280A4, Microbial Sequence Analysis, in the fall of 2022. In this project, I used Illumina and Nanopore sequencing data from a bacterial genome to identify a 

It is written in [Markdown format](https://www.markdownguide.org/basic-syntax/).  

## Step 1: Retrieve Illumina and Nanopore data files

1. Enter server and navigate to directory where data files are stored

2. Copy files into working directory (directory where GitHub repository is located)

```
ADD CP COMMAND
```

## Step 2: Quality check

1. Enter conda environment 

```
conda activate bio_tools
```

2. Check version of fastqc

```
which fastqc
```

3. Check quality of paired-end Illumina data

```
fastqc DATA1.fastq DATA2.fastq
```

4. Check quality of Nanopore data

```
fastqc DATA.fastq
```

## Step 3: Clean Illumina data

1. Check version of cutadapt

```
which cutadapt
```

2. Trim Nextera adapters and low-quality reads

```
cutadapt \
   -a CTGTCTCTTATACACATCT \
   -A CTGTCTCTTATACACATCT \
   -q 30,30 \
   --minimum-length 80 \
   -o DATA1_trimmed.fastq \
   -p DATA2_trimmed.fastq \
   DATA1.fastq \
   DATA2.fastq \
   | tee cutadapt.log
```

## Step 4: Quality recheck

1. Check trimmed Illumina quality

```
fastqc DATA1_trimmed.fastq DATA2_trimmed.fastq
```

2. Compare data before and after trimming

    a. HOW MANY ILLUMINA READS CONTAINED ADAPTERS?

    b. WHAT PERCENT OF READ PAIRS MADE IT THROUGH FILTERING?
    
    c. lINK TO QUALITY REPORT IN GITHUB
    
## Step 5: Assemble Illumina reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_spades_assembly \
   --pe1-1 DATA1_trimmed.fastq \
   --pe1-2 DATA2_trimmed.fastq \
   -m 24 -t 18
```

## Step 6: Assemble Illumina and Nanopore reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_spades_assembly \
   --pe1-1 DATA1_trimmed.fastq \
   --pe1-2 DATA2_trimmed.fastq \
   --nanopore NANOPOREDATA.fastq \
   -m 24 -t 18
```

## Step 7: Compare assemblies

1. Contig lengths
    
    a. WHAT ARE THE LENGTHS OF THE 12 LONGEST CONTIGS OF EACH?
    
2. Contig stats

    a. WHAT IS THE N50 FOR THE 12 LONGEST CONTIGS OF EACH?
    
## Step 8: BLAST contigs
