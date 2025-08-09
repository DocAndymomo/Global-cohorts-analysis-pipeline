input_folder='/mnt/DATA2/TCGA-SKCM/03-fastp/01-data'
output_folder='/mnt/DATA1/GEO/04-bowtie2'

db_bowtie2='/root/1-Biosoft/DB/GRCh38_noalt_as/GRCh38_noalt_as'

# 创建输出文件夹（如果不存在）
mkdir -p $output_folder

# 设置批次大小
batch_size=1

# 计数器
count=0

# 遍历文件夹中的所有fastq文件
for file in "$input_folder"/*.rna_seq.genomic.gdc_realn_fastp_R1.fastq.gz; do
    # 获取文件名和路径
    filename=$(basename "$file")
    # 提取文件名前第一个.前的部分
    file_prefix="${filename%%.*}"
    
    # 指定输出路径和文件名
    output="${file_prefix}_host_remove"
    output_R1="${output}_R1.fastq.gz"
    output_R2="${output}_R2.fastq.gz"
    
    # 检查是否已经存在解析的文件
    if [[ -e "$output_folder/$output_R1" ]] && [[ -e "$output_folder/$output_R2" ]]; then
        echo "跳过文件 ${output_R1} 和 ${output_R2}，因为它们已经存在。"
        continue
    fi
    
    # 构建第二端fastq文件的路径
    input_R2="$input_folder/${file_prefix}.rna_seq.genomic.gdc_realn_fastp_R2.fastq.gz"
    
    # 运行Bowtie2
    bowtie2 -p 45 --quiet -x "$db_bowtie2" -1 "$file" -2 "$input_R2" --un-conc-gz "$output_folder/$output" --no-unal
    
    # 重命名生成的host-sequence free样本
    mv "$output_folder/$output.1" "$output_folder/$output_R1"
    mv "$output_folder/$output.2" "$output_folder/$output_R2"
    
    # 增加计数器
    ((count++))
    
    # 检查计数器是否达到批次大小，如果达到则等待所有进程完成
    if [ "$count" -eq "$batch_size" ]; then
        wait
        count=0
    fi
done

# 确保等待剩余的进程完成
wait
