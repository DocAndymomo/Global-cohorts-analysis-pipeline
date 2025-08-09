#!/bin/bash

# 定义源目录、目标目录和匹配关系文件路径
SOURCE_DIR="/mnt/DATA2/GEO/01-raw/10"
TARGET_DIR="/mnt/DATA2/GEO/01-raw/10cat/"
MATCH_FILE="/mnt/DATA2/GEO/01-raw/10-cat_seq.txt"

# 确保目标目录存在
mkdir -p "$TARGET_DIR"

# 初始化文件描述符列表
declare -A FD_MAP

# 读取匹配文件并处理每一行
while IFS=$'\t' read -r SAMPLE_NAME FASTQ_FILE; do
    # 根据样本名和文件名中的_1或_2后缀，确定目标文件名
    if [[ $FASTQ_FILE =~ _1.fastq.gz ]]; then
        TARGET_FASTQ="$TARGET_DIR/${SAMPLE_NAME}_merged_1.fastq.gz"
    elif [[ $FASTQ_FILE =~ _2.fastq.gz ]]; then
        TARGET_FASTQ="$TARGET_DIR/${SAMPLE_NAME}_merged_2.fastq.gz"
    else
        echo "Warning: File name does not match expecy6ted pattern: $FASTQ_FILE"
        continue
    fi
    
    # 如果目标文件的文件描述符尚未打开，那么打开它
    if [ -z "${FD_MAP[$TARGET_FASTQ]}" ]; then
        exec {FD}>"$TARGET_FASTQ"
        FD_MAP["$TARGET_FASTQ"]=$FD
    else
        FD=${FD_MAP["$TARGET_FASTQ"]}
    fi
    
    # 将文件内容追加到目标文件
    cat "$SOURCE_DIR/$FASTQ_FILE" >&$FD
    
done < "$MATCH_FILE"

# 关闭所有打开的文件描述符
for FD in ${FD_MAP[@]}; do
    exec {FD}>&-
done

echo "All files processed."
