import sys

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 写入输出文件的表头
        outfile.write("SNP\tChromosome\tPosition\ttrait1\n")
        
        # 初始化SNP计数器
        snp_counter = 1
        
        # 处理每一行数据
        for line in infile:
            columns = line.strip().split()
            
            # 确保行有足够的列（至少5列）
            if len(columns) < 5:
                print(f"警告：跳过列数不足的行: {line.strip()}")
                continue
            
            # 提取需要的列
            chromosome = columns[1]  # 第二列（索引1）
            position = columns[2]    # 第三列（索引2）
            trait1 = columns[4]      # 第五列（索引4）
            
            # 写入输出文件
            outfile.write(f"{snp_counter}\t{chromosome}\t{position}\t{trait1}\n")
            
            # SNP计数器递增
            snp_counter += 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("使用方法: python script.py <输入文件> <输出文件>")
        print("示例: python process_data.py input.txt output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_file(input_file, output_file)
    print(f"处理完成！结果已保存到 {output_file}")
