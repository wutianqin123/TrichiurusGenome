#!/bin/bash

# 定义目录变量
input_dir="/mnt/genome10/Lab_Users/wtq_callno/coitree/trim/ten/"
output_dir="/mnt/genome10/Lab_Users/wtq_callno/coitree/getorg/ten/"
reference_path="/mnt/genome10/Lab_Users/wtq_callno/coitree/coipath/ogpath/ten/ten/ten_genomic.fas"

# 确保输出目录存在
mkdir -p "$output_dir"

# 遍历目录中的所有 _R1.fq 文件
for r1_file in "$input_dir"*_R1.fq; do
    # 提取样本名（假设文件名格式为 sample_R1.fq）
    sample_name=$(basename "$r1_file" _R1.fq)
    
    # 构建对应的 _R2.fq 文件路径
    r2_file="$input_dir${sample_name}_R2.fq"
    
    # 检查 _R2.fq 文件是否存在
    if [[ -f "$r2_file" ]]; then
        # 构造输出文件名
        output_file="$output_dir${sample_name}_out"
        
        # 执行命令
        get_organelle_from_reads.py -1 "$r1_file" -2 "$r2_file" -R 10 -s "$reference_path" -k 21,45,65,85,105 -F animal_mt -o "$output_file" -t 5 &
    else
        echo "Error: $r2_file not found for $r1_file"
    fi
done
