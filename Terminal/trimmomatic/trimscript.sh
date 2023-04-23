#!/bin/bash
# Set path to trimmomatic


# Loop through each sample
for sample in SRR10600888 SRR10600877 SRR10600875 SRR10600934 SRR10600933 SRR10600931 SRR10600920 SRR10600919 SRR10600918 SRR10600907 SRR10600906 SRR10600905
do
  # Set input and output file paths
  input1=${sample}_1.fastq.gz
  input2=${sample}_2.fastq.gz
  output1=${sample}_trimmed_1.fastq.gz
  output2=${sample}_trimmed_2.fastq.gz
  
  # Run trimmomatic
  java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 $input1 $input2 $output1 $output2 \
  ILLUMINACLIP: Trimmomatic-0.39/adapters/TrueSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70
done
