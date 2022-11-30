# MIP 280A4 Final Project 2022

## Initial Summary

This report documents the final project I did for MIP 280A4, Microbial Sequence Analysis, in the fall of 2022. In this project, I used Illumina and Nanopore sequencing data from a bacterial genome to identify a 

It is written in [Markdown format](https://www.markdownguide.org/basic-syntax/). 

## Step 1: Clone GitHub repository

1. Clone GitHub repository
```
git clone https://github.com/sophie-kiehl/MIP_280A4_Final_Project.git
```

## Step 2: Retrieve Illumina and Nanopore data files

1. Enter thoth server and navigate to directory where data files are stored

```
ssh skiehl@thoth01.cvmbs.colostate.edu

cd /home/data_for_classes/2022_MIP_280A4/final_project_datasets
```

2. Copy files into directory with GitHub repository and them move into that directory

```
cp Planococcus* ~/MIP_280A4_Final_Project

cd ~/MIP_280A4_Final_Project
```

## Step 3: Quality check

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
fastqc Planococcus_Illumina_R1.fastq Planococcus_Illumina_R2.fastq
```

Quality of first read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R1_fastqc.html

Quality of second read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R2_fastqc.html

4. Check quality of Nanopore data

```
fastqc Planococcus_Nanopore.fastq
```

Quality of read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Nanopore_fastqc.html

5. Create directory to deposit output of fastqc

```
mkdir Quality_Control_and_Trimming

mv ~/MIP_280A4_Final_Project/Planococcus*.html ~/MIP_280A4_Final_Project/Quality_Control_and_Trimming
mv ~/MIP_280A4_Final_Project/Planococcus*.zip ~/MIP_280A4_Final_Project/Quality_Control_and_Trimming
```

6. Push quality reports into GitHub repository

```
cd Quality_Control_and_Trimming

git add Planococcus*.html

git config --global user.email "email address"
git config --global user.name "name"
git config http.postBuffer 524288000

git commit -m "fastqc quality reports for Illumina and Nanopore"

git push origin main
```

7. Push compressed raw NGS data into GitHub repository

```
git add Planococcus*.zip

git commit -m "compressed raw NGS data (Illumina and Nanopore)"

git push origin main
```

## Step 4: Clean Illumina data

1. Check version of cutadapt

```
which cutadapt
```

2. Trim Universal adapters (as this was the most prevalent adapter according to the fastqc quality reports) and remove short (less than 80 bp) and low-quality reads (Phred score beneath 30)

```
cutadapt \
   -a AGATCGGAAGAG \
   -A AGATCGGAAGAG \
   -q 30,30 \
   --minimum-length 80 \
   -o Planococcus_Illumina_R1_trimmed.fastq \
   -p Planococcus_Illumina_R2_trimmed.fastq \
   Planococcus_Illumina_R1.fastq \
   Planococcus_Illumina_R2.fastq \
   | tee ~/MIP_280A4_Final_Project/Quality_Control_and_Trimming/cutadapt.log
```

3. Examine the trimming report

   a. 33.4% of the forward reads and 31.0% of the reverse reads had adapter sequences trimmed
   
   b. 3.0% (25,586) of pairs were remove because they were too short
   
   c. 3.6% of basepairs were removed due to low quality (3,778,458 from the forward read and 11,304,148 from the reverse)
   
   d. In total, 91.1% of base pairs passed filtering

3. Push cutadapt log into GitHub repository

```
cd Quality_Control_and_Trimming

git add cutadapt.log

git commit -m "output of adapter trimming"

git push origin main
```

## Step 5: Quality recheck

1. Check trimmed Illumina quality

```
cd ~/MIP_280A4_Final_Project

fastqc Planococcus_Illumina_R1_trimmed.fastq Planococcus_Illumina_R2_trimmed.fastq
```

2. Push quality report to GitHub repository

```
ADD PUSH COMMAND
```

3. Compare data before and after trimming

    a. HOW MANY ILLUMINA READS CONTAINED ADAPTERS THAT WERE TRIMMED?

    b. WHAT PERCENT OF READ PAIRS MADE IT THROUGH FILTERING?
    
    
## Step 6: Assemble Illumina reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_spades_assembly \
   --pe1-1 DATA1_trimmed.fastq \
   --pe1-2 DATA2_trimmed.fastq \
   -m 24 -t 18
```

2. Push assembly output to GitHub repository

```
ADD PUSH COMMAND
```

## Step 7: Assemble Illumina and Nanopore reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_and_nanopore_spades_assembly \
   --pe1-1 DATA1_trimmed.fastq \
   --pe1-2 DATA2_trimmed.fastq \
   --nanopore NANOPOREDATA.fastq \
   -m 24 -t 18
```

2. Push assembly output to GitHub repository

```
ADD PUSH COMMAND
```

## Step 8: Compare assemblies

1. Number of contigs

    a. HOW MANY CONTIGS ARE IN EACH ASSEMBLY?

2. Number of scaffolds

    a. HOW MANY SCAFFOLDS ARE IN EACH ASSEMBLY?

3. Contig lengths
    
    a. WHAT ARE THE LENGTHS OF THE 12 LONGEST CONTIGS OF EACH?
    
4. Contig stats

    a. WHAT IS THE N50 FOR THE 12 LONGEST CONTIGS OF EACH?
    
## Step 9: BLAST contigs

1. If there are a larger number of contigs, extract top 12 contigs from Illumina and Nanopore assembly

```
ADD EXTRACT COMMAND
```

2. Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

3. Record search parameters

4. Record top 5 hits and percent identity for each of the first 12 contigs

## Step 10: Genome alignment

