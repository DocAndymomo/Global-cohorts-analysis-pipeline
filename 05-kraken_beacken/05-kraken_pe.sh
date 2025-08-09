#!/bin/bash
#以上代码可以合并成： 
input_folder="/mnt/DATA1/TCGA/04-bowtie2" 
output_folder="/mnt/DATA1/TCGA/05-ku" 
report_folder=$output_folder 
db_ku="/home/rstudio/SAHMI/k2_pluspf_20230605" 

# 创建输出文件夹和报告文件夹
mkdir -p $output_folder mkdir -p $report_folder
# 获取input_folder下所有fastq_gz文件的列表
files=$(find $input_folder -name "*_host_remove_R1.fastq.gz") 
count=0 
total=$(echo $files | wc -w)
# 迭代处理每个文件的R1和R
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
    output_file_name=${file_base//"host_remove"/}
    # 构建输出文件路径和报告文件路径
    output_file=$output_folder/${output_file_name}.kraken.output.txt 
    report=$report_folder/${output_file_name}.kraken.report.txt 
    report_std=$report_folder/${output_file_name}.kraken.report.std.txt 
    mpa_output=$report_folder/${output_file_name}.kraken.mpa.txt
    # 在处理之前，检查结果文件是否存在
    if [ -f $report_std ]; then 
        echo "$report_std 已存在，跳过处理。" 
      continue # 跳过处理该样本 
    fi
    
    # 运行kraken2命令
    /home/rstudio/SAHMI/kraken2-2.1.3/kraken2 --db $db_ku --threads 32 --report-minimizer-data --report $report --use-names --output $output_file --paired $input_fq1 $input_fq2
    
    cut -f1-3,6-8 $report > $report_std
    
    python "${code_path}/kreport2mpa.py" -r $report_std -o $mpa_output
done

