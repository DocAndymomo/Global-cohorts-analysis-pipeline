input_folder='/mnt/DATA1/TCGA/03-fastp/01-data'
output_folder='/mnt/DATA1/TCGA/04-bowtie2sen-38/'
db_bowtie2='/root/1-Biosoft/DB/GRCh38_noalt_as/GRCh38_noalt_as'

# 创建输出文件夹（如果不存在）
mkdir -p "$output_folder"


# 遍历文件夹中的所有fastq文件
for file in  `ls -r "$input_folder"/*_R1.fastq.gz`; do
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
    
    # 运行Bowtie2，并将输出重定向到 /dev/null
    echo "Running bowtie2 for $file and $input_R2..."
    /root/1-Biosoft/bowtie2-2.5.3-linux-x86_64/bowtie2 -p 40 --quiet -x "$db_bowtie2" -1 "$file" -2 "$input_R2" --un-conc-gz "$output_folder/$output" --no-unal --end-to-end --very-sensitive -k 16 --np 1 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05 > /dev/null 2>&1 &
    
    # 等待 bowtie2 任务完成
    wait
    
    # 检查并重命名生成的host-sequence free样本
    if [[ -e "$output_folder/${output}.1" ]] && [[ -e "$output_folder/${output}.2" ]]; then
        mv "$output_folder/${output}.1" "$output_folder/$output_R1" > /dev/null 2>&1
        mv "$output_folder/${output}.2" "$output_folder/$output_R2" > /dev/null 2>&1
        echo "Processed and moved files for $file_prefix"
    else
        echo "Error: Bowtie2 did not generate expected output files for $file_prefix"
    fi
    
done

# 确保等待剩余的进程完成
wait
