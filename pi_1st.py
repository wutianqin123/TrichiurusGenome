import sys
import os

def process_file(input_file):
    # 定义四个输出文件
    output_files = {
        'popA': open('popA_output.txt', 'w'),
        'popB': open('popB_output.txt', 'w'),
        'popC': open('popC_output.txt', 'w'),
        'popD': open('popD_output.txt', 'w')
    }
    
    # 读取并处理数据
    data = {'popA': [], 'popB': [], 'popC': [], 'popD': []}
    
    with open(input_file, 'r') as f:
        for line in f:
            # 跳过空行
            if not line.strip():
                continue
                
            parts = line.strip().split()
            if len(parts) < 5:
                continue
                
            pop = parts[0]
            chromosome = parts[1]
            position = parts[2]
            # 处理第五列，将NA替换为0
            value = '0' if parts[4].upper() == 'NA' else parts[4]
            
            # 只处理染色体1-24
            if chromosome.isdigit() and 1 <= int(chromosome) <= 24:
                # 保存数据用于排序
                data[pop].append((int(chromosome), int(position), value, line.strip()))
    
    # 对每个群体的数据进行排序并写入文件
    for pop, records in data.items():
        # 按染色体(第一优先级)和位点(第二优先级)排序
        records.sort(key=lambda x: (x[0], x[1]))
        
        # 写入排序后的数据
        for chrom, pos, val, original_line in records:
            # 分割原始行并替换第五列
            parts = original_line.split()
            parts[4] = val  # 替换处理后的值
            output_files[pop].write('\t'.join(parts) + '\n')
    
    # 关闭所有文件
    for file in output_files.values():
        file.close()
    
    print("处理完成！文件已保存为:")
    for pop in ['popA', 'popB', 'popC', 'popD']:
        print(f"- {pop}_output.txt")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使用方法: python script.py <输入文件名>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"错误: 文件 '{input_file}' 不存在")
        sys.exit(1)
    
    process_file(input_file)
