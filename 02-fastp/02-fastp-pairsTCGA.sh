#!/bin/bash
input_dir="/mnt/usbshare1/02-melanoma_GEO_backup"
output_dir="/mnt/DATA1/01-Healthy/03-fastp/01-data/"  # 添加尾部斜杠确保目录结构正确
html_dir="/mnt/DATA1/01-Healthy/03-fastp/02-html/"
mkdir -p "$output_dir"
count=0 # 当前正在运行的进程数

for forward_file in "$input_dir"/*_1.fastq.gz; do
  # 提取文件名和扩展名
  filename=$(basename "$forward_file")
  filename_no_ext="${filename%_1.fastq.gz}"
  # 生成反向文件名
  reverse_file="${input_dir}/${filename_no_ext}_2.fastq.gz"
  # 生成输出文件路径，文件名中加入fastp
  output_forward="${output_dir}${filename_no_ext}_fastp_1.fastq.gz"
  output_reverse="${output_dir}${filename_no_ext}_fastp_2.fastq.gz"
  
  # 检查输出文件是否已存在，如果存在则跳过循环
  if [ -f "$output_forward" ] && [ -f "$output_reverse" ]; then
    echo "Skipped: $forward_file"
    echo "Skipped: $reverse_file"
    continue
  fi
  
  # 执行 fastp 命令，并放到后台
  fastp -i "$forward_file" -o "$output_forward" -I "$reverse_file" -O "$output_reverse" -Q --thread=48 --length_required=45 –stdout –interleaved --compression=6 --html "${html_dir}${filename_no_ext}_fastp.html" &
  echo "Processed: $forward_file"
  echo "Processed: $reverse_file"
  ((count++)) # 当前正在运行的进程数加1
  
  if [ $count -ge 80 ]; then
    wait # 等待前面的进程完成
    count=0 # 重置计数器
  fi

done

wait  # 等待剩余的进程完成