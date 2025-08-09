#!/bin/bash

# 定义输入、输出和报告文件夹
input_folder="/mnt/DATA1/GEO-1/04-hr2"
output_folder="/mnt/DATA1/GEO-1/05-ku"
report_folder=$output_folder
db_ku="/mnt/DATA4/DB_3/k2_database/01-Microbial2023"

# 创建输出文件夹和报告文件夹
mkdir -p $output_folder
mkdir -p $report_folder

# 获取input_folder下所有fastq_gz文件的列表
files=$(find $input_folder -name "*_chm_remove_R1.fastq.gz")
count=0
total=$(echo $files | wc -w)

# 运行并发任务数
max_jobs=1
current_jobs=0

# 迭代处理每个文件的R1和R2
for file in $files; do
    # 提取文件名和路径
    count=$((count+1))
    echo "Processing file $count of $total: $file"
    file_name=$(basename $file)
    file_path=$(dirname $file)
    
    # 构建R2文件路径
    file_base="${file_name%_R1.fastq.gz}"
    input_fq1=$file_path/$file_base"_R1.fastq.gz"
    input_fq2=$file_path/$file_base"_R2.fastq.gz"
    
    # 生成不包含"host_remove"字符的输出文件名
    output_file_name=${file_base//"chm_remove"/}
    
    # 构建输出文件路径和报告文件路径
    report=$report_folder/${output_file_name}.kraken.report.txt
    report_std=$report_folder/${output_file_name}.kraken.report.std.txt
    mpa_output=$report_folder/${output_file_name}.kraken.mpa.txt
    
    # 在处理之前，检查结果文件是否存在
    if [ -f $report_std ]; then
        echo "$report_std 已存在，跳过处理。"
        continue # 跳过处理该样本
    fi
    
    # 运行krakenuniq命令
    krakenuniq \
        --db $db_ku \
        --threads 40 \
        --preload \
        --report-file $report \
        --paired $input_fq1 $input_fq2 &
    
    current_jobs=$((current_jobs+1))
    
    # 等待当前任务完成
    if [ $current_jobs -ge $max_jobs ]; then
        wait
        current_jobs=0
    fi

    cut -f1-3,6-8 $report > $report_std
    
    python "/root/miniconda3/bin/kreport2mpa.py" -r $report_std -o $mpa_output
done

# 等待所有后台任务完成
wait
