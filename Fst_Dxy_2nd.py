import sys

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 读取表头行
        header = infile.readline().strip().split()
        
        # 检查输入文件是否有足够的列
        if len(header) < 6:
            print("错误：输入文件至少需要6列")
            return
        
        # 写入输出文件的表头
        outfile.write("SNP\tChromosome\tPosition\ttrait1\n")
        
        # 初始化SNP计数器
        snp_counter = 1
        
        # 处理每一行数据
        for line in infile:
            columns = line.strip().split()
            
            # 确保行有足够的列
            if len(columns) < 6:
                continue
            
            # 提取需要的列
            chromosome = columns[2]  # 第三列
            position = columns[3]    # 第四列
            trait1 = columns[5]      # 第六列
            
            # 写入输出文件
            outfile.write(f"{snp_counter}\t{chromosome}\t{position}\t{trait1}\n")
            
            # SNP计数器递增
            snp_counter += 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("使用方法: python script.py <输入文件> <输出文件>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_file(input_file, output_file)
    print(f"处理完成！结果已保存到 {output_file}")
