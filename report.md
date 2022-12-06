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
fastqc -t 24 Planococcus_Nanopore.fastq.gz
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
   --nanopore Planococcus_Nanopore.fastq.gz \
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

3. Find the length of the 12 longest contigs and scaffolds for each assembly

```
grep NODE ./paired_spades_assembly/contigs.fasta | head -12
grep NODE ./paired_and_nanopore_spades_assembly/contigs.fasta | head -12

grep NODE ./paired_spades_assembly/scaffolds.fasta | head -12
grep NODE ./paired_and_nanopore_spades_assembly/scaffolds.fasta | head -12
```
The length of each is in the header of each node

| | Illumina only contigs | Illumina and Nanopore contigs | Illumina only scaffolds | Illumina and Nanopore scaffolds |
| --- | --- | --- | --- | --- |
| Number | 32 | 12 | 30 | 11 |
| N50 | 453155 | 1073186 | 453155 | 1073186 |
| 1 length | 1662249 | 1447436 | 1662249 | 1793883 |
| 2 length | 453155 | 1073186 | 453155 | 1447436 |
| 3 length | 437964 | 718909 | 437964 | 234788 |
| 4 length | 346434 | 234788 | 386565 | 48414 |
| 5 length | 174533 | 48414 | 174533 | 12386 |
| 6 length | 98084 | 12386 | 98084 | 7491 |
| 7 length | 90319 | 7491 | 90319 | 400 |
| 8 length | 65255 | 400 | 65255 | 399 |
| 9 length | 51313 | 399 | 51313 | 383 |
| 10 length | 48414 | 383 | 48414 | 186 |
| 11 length | 40031 | 186 | 15155 | 128 |
| 12 length | 15155 | 128 | 12386 | N/A |

    
## Step 9: BLAST contigs
1. Extract the first twelve scaffolds from the Illumina and Nanopore assembly

```
cd paired_and_nanopore_spades_assembly

seqtk seq -A scaffolds.fasta | head -24 > first_12_scaffolds.fasta
```

2. Retrive the first_12_scaffolds.fasta file from the server

Open a new terminal window and enter the directory you would like to place the file in

```
sftp "thoth username"

cd MIP_280A4_Final_Project/paired_and_nanopore_spades_assembly

get first_12_scaffolds_fasta
```

3. Open first_12_contigs.fasta in a plain text editor (such as BBEdit)

4. Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

5. Record search parameters

6. Record top 5 hits and percent identity for each of the first 12 contigs

| Scaffold (run) | 1 | % | 2 | % | 3 | % | 4 | % | 5 | % |
| --- | ---| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 (1) Description Graphics | Planococcus sp. chr | 99 | Planococcus rifietoensis chr | 99 | Planococcus plakortidis chr | 99 | Planococcus maritimus chr | 99 | Planococcus maritimus chr | 99 |
| 1 (2) Description Graphics | Planococcus sp. chr | 93.05 | Planococcus rifietoensis chr | 92.51 | Planococcus maritimus | 91.63 | Planococcus plakortidis | 88.90 | Planococcus maritimus | 91.73 |
| 1 (3) Description Graphics | Planococcus rifietoensis chr | 96.79 | Planococcus plakortidis chr | 92.38 | Planococcus sp. chr | 95.20 | Planococcus maritimus chr | 88.68 | Planococcus kocurii chr | 80.0 |
| 2 (1) Description Graphics | Planococcus rifietoensis chr | 94.86 | Planococcus sp. | 89.28 | Planococcus maritimus chr | 90.76 | Planococcus plakortidis chr | 91.50 | Planococcus glaciei chr | 76.68 |
| 2 (2) Description Graphics | Planococcus sp. chr | 92.04 | Planococcus rifietoensis chr | 93.41 | Planococcus maritimus chr | 89.68 | Planococcus plakortidis chr | 88.94 | Planococcus maritimus chr | 86.04 |
| 2 (3) Description Graphics | Planococcus rifietoensis chr | 88.72 | Planococcus plakortidis chr| 96.03 | Planococcus sp. chr | 91.29 | Planococcus maritimus chr | 86.11 | Planococcus maritmus chr | 80.20 |
| 3 (1) Description Graphics | Planococcus rifietoensis chr | 90.37 | Planococcus sp. chr | 90.39 | Planococcus plakortidis chr | 89.07 | Planococcus maritimus chr | 85.77 | Planococcus maritimus chr | 84.55 |
| 3 (2) Description Graphics | ? | ? | ? | ? | ? | ? | ? | ? | ? | ? |
| 3 (3) Description Graphics | ? | ? | ? | ? | ? | ? | ? | ? | ? | ? |
| 4 Description Graphics | Planococcus maritimus plasmid | 90.35 | Planococcus sp. plasmid | 92.25 | Planococcus faecalis chr | 90.86 | Planococcus maritimus chr | 94.17 | Planococcus antarcticus | 92.90 |
| 5 [Description](https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/BLAST_results/5_Description.png) [Graphics](https://github.com/sophie-kiehl/MIP_280A4_Final_Project/blob/main/BLAST_results/5_Graphics.png) | Planococcus citreus plasmid | 95.35 | Planococcus sp. plasmid  | 94.69 | Planococcus sp. plasmid | 94.60 | Staphylococcus sciuri plasmid | 89.80 | Mammaliicoccus sciuri plasmid | 89.80 |
| 6 Description Graphics | Planococcus citreus plasmid | 94.69 | Planococcus sp. plasmid | 94.21 | Planococcus sp. plasmid | 93.61 | Staphylococcus scirui plasmid | 94.09 | Planococcus kocurii plasmid | 87.33 |

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
