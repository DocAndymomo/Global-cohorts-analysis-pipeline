#!/bin/bash

input_folder="/mnt/DATA1/GEO/01-data/1/"
output_folder="/mnt/DATA2/GEO/02-fastq/1-2479/"
temp_folder="/mnt/DATA2/temp/"

# 创建输出文件夹（如果不存在）
mkdir -p "$output_folder"
# 创建临时文件夹（如果不存在）
mkdir -p "$temp_folder"

# 查找所有 .sra 文件，包括子目录中的文件
find "$input_folder" -name "*.sra" | while read -r sra_file; do
    # 检查是否是文件
    if [[ -f "$sra_file" ]]; then
        # 获取文件名（不含扩展名）
        filename=$(basename -- "$sra_file")
        filename_noext="${filename%.sra}"

        # 指定未压缩的输出文件路径
        output_file_prefix="$output_folder/${filename_noext}"

        # 检查输出文件是否存在
        output_exists=false
        for fastq_file in "${output_file_prefix}"*.fastq.gz; do
            if [[ -f "$fastq_file" ]]; then
                output_exists=true
                break
            fi
        done

        if [ "$output_exists" = true ]; then
            echo "Skipped: $sra_file (output files already exist)"
            continue
        fi

        # 执行转换，并指定输出目录和临时文件目录
        /root/1-Biosoft/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump "$sra_file" --split-3 -p -e 96 -O "$output_folder" --temp "$temp_folder"

        # 压缩输出文件。假设 `--split-3` 会产生 1 到 3 个 fastq 文件
        for fastq_file in "${output_file_prefix}"*.fastq; do
            if [[ -f "$fastq_file" ]]; then
                pigz -p 96 "$fastq_file"
            fi
        done

        echo "Processed: $sra_file"
    fi
done

# 统计未完成转换的 .sra 文件
echo "未完成转换的 .sra 文件："
find "$input_folder" -name "*.sra" | while read -r sra_file; do
    filename=$(basename -- "$sra_file")
    filename_noext="${filename%.sra}"
    output_file_prefix="$output_folder/${filename_noext}"
    
    # 检查是否存在对应的 .fastq.gz 文件
    fastq_files_count=$(find "${output_file_prefix}"*.fastq.gz 2>/dev/null | wc -l)
    if [ "$fastq_files_count" -eq 0 ]; then
        echo "$sra_file"
    fi
done

