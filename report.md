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

There were 844,025 forward reads in total. There was slightly abnormal distribution of GC content and, as expected, some of the reads included adapters, specifically the Illumina Universal adapter.

<img width="859" alt="Warning on GC content distribution" src="https://user-images.githubusercontent.com/115187825/204867185-162de0b7-be05-4682-997c-9ac33bda156c.png">

<img width="901" alt="Warning on adapter content" src="https://user-images.githubusercontent.com/115187825/204867226-96286b9e-dc13-4ff2-9071-8546ec001c7a.png">

Quality of second read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R2_fastqc.html

There were 844,025 reverse reads in total. In addition to an abnormal GC content distribution and the presence of adapters, the reverse reads also showed warnings on per base sequence quality and on per base sequence content.

<img width="930" alt="Warning on per base sequence quality" src="https://user-images.githubusercontent.com/115187825/204868491-bb77b74f-e0ea-4031-853b-4e8530cc48ed.png">

<img width="929" alt="Warning on per base sequence content" src="https://user-images.githubusercontent.com/115187825/204868512-1d4b0f77-ca93-46c1-9c40-ec84edb06377.png">

4. Check quality of Nanopore data

```
fastqc Planococcus_Nanopore.fastq
```

Quality of read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Nanopore_fastqc.html

There were 1371 reads. The quality report showed errors on per base sequence quality, per base sequence content, and per sequence GC content and warnings on per sequence quality scores and sequence length distribution. This is not unexepcted quality for Nanopore reads.

<img width="872" alt="Error on per base sequence quality" src="https://user-images.githubusercontent.com/115187825/204870502-4722776d-fd0f-4b58-bd48-2d83244f3881.png">
on.

<img width="871" alt="Error on per base sequence content" src="https://user-images.githubusercontent.com/115187825/204870613-8cea2ba2-f858-4122-8ff0-16ba00cf80f8.png">


5. Create directory to deposit output of fastqc

```
mkdir Quality_Control_and_Trimming

mv Planococcus*.html ./Quality_Control_and_Trimming
mv Planococcus*.zip ./Quality_Control_and_Trimming
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

4. Push cutadapt log into GitHub repository

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

2. Move quality report Quality Control direcotyr and push to GitHub repository

```
mv Planococcus*trimmed_fastqc.html ./Quality_Control_and_Trimming
mv Plnaococcus*trimmed_fastqc.zip ./Quality_Control_and_Trimming

cd Quality_Control_and_Trimming

git add Planococcus*trimmed_fastqc.html
git add Planococcus*trimmed_fastqc.zip

git commit -m "quality reports after trimming and filtering"

git push origin main
```

3. Compare data before and after trimming

Quality of first read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R1_trimmed_fastqc.html

There were 818,438 reads in total. All of the adapter content was successfully removed. The abnormal GC content distribution remained. Additionally, this quality report showed a warning on sequence length distribution that was not present on the report pre-filtering.

<img width="853" alt="Warning on sequence lenght distribution" src="https://user-images.githubusercontent.com/115187825/204867695-9fc99918-1636-4569-a069-8fd1b5fe855f.png">

Quality of second read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R2_trimmed_fastqc.html

There were 818,439 reverse reads after trimming. Otherwise, the results from this were the same as the results for trimming the R1 reads.
    
    
## Step 6: Assemble Illumina reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_spades_assembly \
   --pe1-1 Planococcus_Illumina_R1_trimmed.fastq \
   --pe1-2 Planococcus_Illumina_R2_trimmed.fastq \
   -m 24 -t 18
```

2. Extract the 12 largest contigs

```
seqtk seq -A contigs.fasta | head -24 > first_12_contigs.fasta
```

2. Push assembly output to GitHub repository

```
ADD PUSH COMMAND
```

## Step 7: Assemble Illumina and Nanopore reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_and_nanopore_spades_assembly \
   --pe1-1 Planococcus_Illumina_R1_trimmed.fastq \
   --pe1-2 Planococcus_Illumina_R2_trimmed.fastq \
   --nanopore Planococcus_Nanopore.fastq \
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

