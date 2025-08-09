#!/bin/bash

# 输入文件夹
input_dir="/media/h11/1CDE-C172/TCGA-SKCM_genomic"
# 输出文件夹
output_dir="/media/h11/DATA2/TCGA-SKCM/02-fastq"

# 确保输出文件夹存在
mkdir -p "$output_dir"

# 初始化计数器
count=0
max_jobs=4

# 遍历输入文件夹中的所有BAM文件
for bam_file in "$input_dir"/*.bam; do
  # 获取不带路径的文件名
  filename=$(basename -- "$bam_file")
  # 移除扩展名，准备输出文件名
  base_name="${filename%.bam}"
  
  # 定义输出的 FASTQ 文件路径
  fastq_r1="$output_dir/${base_name}_R1.fastq.gz"
  fastq_r2="$output_dir/${base_name}_R2.fastq.gz"

  # 检查输出文件是否已经存在
  if [[ -f "$fastq_r1" && -f "$fastq_r2" ]]; then
    echo "Output files for $filename already exist, skipping conversion."
    continue
  fi

  # 使用samtools进行排序，并通过管道直接传递给samtools fastq
  (
    samtools sort -n -@ 32 "$bam_file" | \
    samtools fastq -@ 48 - \
      -1 "$fastq_r1" \
      -2 "$fastq_r2" \
      -0 /dev/null -s /dev/null -n
  ) &

  # 计数器递增
  ((count++))

  # 如果达到了最大并行任务数，等待所有后台任务完成
  if ((count == max_jobs)); then
    wait
    count=0
  fi
done

# 等待所有后台任务完成
wait

echo "Conversion completed."
