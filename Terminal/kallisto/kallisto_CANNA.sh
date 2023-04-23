#!/bin/bash
#This script peforms allignment of the Cannabis Sativa Female cDNA transcriptome with Kallisto first, we build an index from our reference fasta file 
#The reference plant transcriptome files are sourced from here: https://plants.ensembl.org/info/data/ftp/index.html 

kallisto index -i CANNA Cannabis_sativa_female.cs10.cdna.all.fa

#now we will map the reads to the indexed reference transcriptome
#we will also create log files to save the output shown in the terminal
#flower & Trichomes beginning from  00920
./kallisto/kallisto quant -i CANNA -o SRR10600888 SRR10600888_trimmed_1.fastq.gz SRR10600888_trimmed_2.fastq.gz &> kallisto_CANNA_888.log
./kallisto/kallisto quant -i CANNA -o SRR10600877 SRR10600877_trimmed_1.fastq.gz SRR10600877_trimmed_2.fastq.gz &> kallisto_CANNA_877.log
./kallisto/kallisto quant -i CANNA -o SRR10600875 SRR10600875_trimmed_1.fastq.gz SRR10600875_trimmed_2.fastq.gz &> kallisto_CANNA_875.log
./kallisto/kallisto quant -i CANNA -o SRR10600934 SRR10600934_trimmed_1.fastq.gz SRR10600934_trimmed_2.fastq.gz &> kallisto_CANNA_934.log
./kallisto/kallisto quant -i CANNA -o SRR10600933 SRR10600933_trimmed_1.fastq.gz SRR10600933_trimmed_2.fastq.gz &> kallisto_CANNA_933.log
./kallisto/kallisto quant -i CANNA -o SRR10600931 SRR10600931_trimmed_1.fastq.gz SRR10600931_trimmed_2.fastq.gz &> kallisto_CANNA_931.log
./kallisto/kallisto quant -i CANNA -o SRR10600920 SRR10600920_trimmed_1.fastq.gz SRR10600920_trimmed_2.fastq.gz &> kallisto_CANNA_920.log
./kallisto/kallisto quant -i CANNA -o SRR10600919 SRR10600919_trimmed_1.fastq.gz SRR10600919_trimmed_2.fastq.gz &> kallisto_CANNA_919.log
./kallisto/kallisto quant -i CANNA -o SRR10600918 SRR10600918_trimmed_1.fastq.gz SRR10600918_trimmed_2.fastq.gz &> kallisto_CANNA_918.log
./kallisto/kallisto quant -i CANNA -o SRR10600907 SRR10600907_trimmed_1.fastq.gz SRR10600907_trimmed_2.fastq.gz &> kallisto_CANNA_907.log
./kallisto/kallisto quant -i CANNA -o SRR10600906 SRR10600906_trimmed_1.fastq.gz SRR10600906_trimmed_2.fastq.gz &> kallisto_CANNA_906.log
./kallisto/kallisto quant -i CANNA -o SRR10600905 SRR10600905_trimmed_1.fastq.gz SRR10600905_trimmed_2.fastq.gz &> kallisto_CANNA_905.log


