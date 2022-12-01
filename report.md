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

2. Check quality of paired-end Illumina data
```
fastqc Planococcus_Illumina_R1.fastq Planococcus_Illumina_R2.fastq
```

Quality of first read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Nanopore_fastqc.html

There were 844,025 forward reads in total. There was slightly abnormal distribution of GC content and, as expected, some of the reads included adapters, specifically the Illumina Universal adapter.

<img width="859" alt="Warning on GC content distribution" src="https://user-images.githubusercontent.com/115187825/204867185-162de0b7-be05-4682-997c-9ac33bda156c.png">

<img width="901" alt="Warning on adapter content" src="https://user-images.githubusercontent.com/115187825/204867226-96286b9e-dc13-4ff2-9071-8546ec001c7a.png">

Quality of second read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Illumina_R2_fastqc.html

There were 844,025 reverse reads in total. In addition to an abnormal GC content distribution and the presence of adapters, the reverse reads also showed warnings on per base sequence quality and on per base sequence content.

<img width="930" alt="Warning on per base sequence quality" src="https://user-images.githubusercontent.com/115187825/204868491-bb77b74f-e0ea-4031-853b-4e8530cc48ed.png">

<img width="929" alt="Warning on per base sequence content" src="https://user-images.githubusercontent.com/115187825/204868512-1d4b0f77-ca93-46c1-9c40-ec84edb06377.png">

3. Check quality of Nanopore data

```
fastqc -t 24 Planococcus_Nanopore.fastq
```

Quality of read: https://htmlpreview.github.io/?https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/Quality_Control_and_Trimming/Planococcus_Nanopore_fastqc.html

There were 6045 reads. The quality report showed errors on per base sequence quality, per sequence quality scores, per base sequence content, and per sequence GC content and a warning on sequence length distribution. This is not unexepcted quality for Nanopore reads.

<img width="1011" alt="Screen Shot 2022-12-01 at 10 30 25 AM" src="https://user-images.githubusercontent.com/115187825/205124505-2597736f-c0d3-4a7c-bf20-bfc8527236cc.png">

<img width="1015" alt="Screen Shot 2022-12-01 at 10 30 37 AM" src="https://user-images.githubusercontent.com/115187825/205124555-ef5b07ec-5ffc-4240-8c79-16f3ded3981c.png">


4. Create directory to deposit output of fastqc

```
mkdir Quality_Control_and_Trimming

mv Planococcus*.html ./Quality_Control_and_Trimming
mv Planococcus*.zip ./Quality_Control_and_Trimming
```

5. Push quality reports into GitHub repository

```
cd Quality_Control_and_Trimming

git add Planococcus*.html

git config --global user.email "email address"
git config --global user.name "name"
git config http.postBuffer 524288000

git commit -m "fastqc quality reports for Illumina and Nanopore"

git push origin main
```

6. Push compressed raw NGS data into GitHub repository

```
git add Planococcus*.zip

git commit -m "compressed raw NGS data (Illumina and Nanopore)"

git push origin main
```

## Step 4: Clean Illumina data

1. Trim Universal adapters (as this was the most prevalent adapter according to the fastqc quality reports) and remove short (less than 80 bp) and low-quality reads (Phred score beneath 30)

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

2. Examine the trimming report

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

## Step 7: Assemble Illumina and Nanopore reads

1. Create de novo assembly using SPAdes assembler (in conda environment)

```
spades.py   -o paired_and_nanopore_spades_assembly \
   --pe1-1 Planococcus_Illumina_R1_trimmed.fastq \
   --pe1-2 Planococcus_Illumina_R2_trimmed.fastq \
   --nanopore Planococcus_Nanopore.fastq \
   -m 24 -t 18
```

## Step 8: Compare assemblies

1. Find the number of contigs and scaffolds for each assembly

```
grep -c NODE ./paired_spades_assembly/contigs.fasta
grep -c NODE ./paired_and_nanopore_spades_assembly/contigs.fasta

grep -c NODE ./paired_spades_assembly/scaffolds.fasta
grep -c NODE ./paired_and_nanopore_spades_assembly/scaffolds.fasta
```

2. Find the N50 values of the contigs and scaffolds for each assembly
```
conda install -c bioconda quast

quast.py ./paired_spades_assembly/contigs.fasta \
        -o Illumina_contigs_quast_test_output
        
quast.py ./paired_and_nanopore_spades_assembly/contigs.fasta \
        -o Illumina_and_nanopore_contigs_quast_test_output
        
quast.py ./paired_spades_assembly/scaffolds.fasta \
        -o Illumina_scaffolds_quast_test_output
        
quast.py ./paired_and_nanopore_spades_assembly/contigs.fasta \
        -o Illumina_and_nanopore_scaffolds_quast_test_output
```

3. Find the length of the 12 longest contigs for each assembly

```
grep NODE ./paired_spades_assembly/contigs.fasta | head -12
grep NODE ./paired_and_nanopore_spades_assembly/contigs.fasta | head -12

grep NODE ./paired_spades_assembly/scaffolds.fasta | head -12
grep NODE ./paired_and_nanopore_spades_assembly/scaffolds.fasta | head -12
```
The length of each contig is in the header of each node

| | Illumina only contigs | Illumina and Nanopore contigs | Illumina only scaffolds | Illumina and Nanopore scaffolds |
| --- | --- | --- | --- | --- |
| Number | 32 | 24 | 30 | 23 |
| N50 | 453155 | 716081 | 453155 | 716081 |
| 1 length | 1662249 | 1662206 | 1662249 | 1662206 |
| 2 length | 453155 | 716081 | 453155 | 716081 |
| 3 length | 437964 | 454644 | 437964 | 548070 |
| 4 length | 346434 | 388083 | 386565 | 388083 |
| 5 length | 174533 | 93258 | 174533 | 65255 |
| 6 length | 98084 | 65255 | 98084 | 51313 |
| 7 length | 90319 | 51313 | 90319 | 48414 |
| 8 length | 65255 | 48414 | 65255 | 15214 |
| 9 length | 51313 | 15214 | 51313 | 12386 |
| 10 length | 48414 | 12386 | 48414 | 7491 |
| 11 length | 40031 | 7491 | 15155 | 429 |
| 12 length | 15155 | 429 | 12386 | 419 |

    
## Step 9: BLAST contigs
4. Extract the first contig

```
head -n 27706 contigs.fasta > contig_1.fasta
```

1. If there are a larger number of contigs, extract top 12 contigs from Illumina and Nanopore assembly

```
ADD EXTRACT COMMAND
```

2. Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

3. Record search parameters

4. Record top 5 hits and percent identity for each of the first 12 contigs

## Step 10: Genome alignment

1. Build an index

```
bowtie2-build ./paired_and_nanopore_spades_assembly/scaffolds.fasta Illumina_and_nanopore_index
```

2. Align the reads to the scaffolds

```
bowtie2 -x Illumina_and_nanopore_index \
   -1 Planococcus_Illumina_R1_trimmed.fastq \
   -2 Planococcus_Illumina_R2_trimmed.fastq \
   --no-unal \
   --threads 8 \
   -S Illumina_reads_mapped_to_scaffolds.sam
```

3. Convert the sam file to a bam file

```
samtools view -b Illumina_reads_mapped_to_scaffolds.sam > Illumina_reads_mapped_to_scaffolds.bam
```
4. Sort the bam files

```
samtools sort -T tmp -O 'bam' Illumina_reads_mapped_to_scaffolds.bam  > Illumina_reads_mapped_to_scaffolds.sorted.bam
```

5. Use samtools to find the average coverage of the Illumina reads to the Illumina and Nanopore scaffolds

```
samtools depth Illumina_reads_mapped_to_scaffolds.sorted.bam > Illumina_reads_mapped_to_scaffolds_coverage.txt
```
