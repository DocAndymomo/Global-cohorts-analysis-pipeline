#!/bin/bash

input_folder="/mnt/DATA2/GEO/03-fastp/test1"
output_folder="/mnt/DATA2/GEO/03-fastp/test2"
mkdir -p "$output_folder"

# 遍历所有 .fastq.gz 文件
find "$input_folder" -type f -name "*.fastq.gz" | while read fastq_file; do
    # 运行 FastQC
    fastqc "$fastq_file" -o "$output_folder"
    echo "Processed: $fastq_file"
done
