import pandas as pd
import numpy as np
import os

def process_fst_file(input_file):
    """
    处理FST文件，将第六列的NA和负数值替换为0
    """
    # 读取文件
    df = pd.read_csv(input_file, sep='\t')
    
    # 检查是否有第六列
    if len(df.columns) < 6:
        print("警告: 文件列数不足6列，无法处理第六列")
    else:
        # 处理第六列的NA和负数值
        sixth_column = df.columns[5]  # 第六列的名称
        
        # 统计处理前的NA和负数值数量
        na_count_before = df[sixth_column].isna().sum()
        negative_count_before = (pd.to_numeric(df[sixth_column], errors='coerce') < 0).sum()
        
        print(f"第六列 '{sixth_column}' 处理前: {na_count_before} 个NA值, {negative_count_before} 个负数值")
        
        # 将NA替换为0
        df[sixth_column] = df[sixth_column].fillna(0)
        
        # 将字符串"NA"替换为0（如果存在）
        if df[sixth_column].dtype == 'object':
            df[sixth_column] = df[sixth_column].replace('NA', 0)
            df[sixth_column] = df[sixth_column].replace('na', 0)
            df[sixth_column] = df[sixth_column].replace('NaN', 0)
            df[sixth_column] = df[sixth_column].replace('nan', 0)
        
        # 转换为数值类型，处理可能的其他非数值
        df[sixth_column] = pd.to_numeric(df[sixth_column], errors='coerce')
        
        # 将负数值替换为0
        df[sixth_column] = np.where(df[sixth_column] < 0, 0, df[sixth_column])
        
        # 再次填充可能因转换产生的NA值
        df[sixth_column] = df[sixth_column].fillna(0)
        
        # 统计处理后的情况
        na_count_after = df[sixth_column].isna().sum()
        negative_count_after = (df[sixth_column] < 0).sum()
        
        print(f"第六列 '{sixth_column}' 处理后: {na_count_after} 个NA值, {negative_count_after} 个负数值")
    
    # 定义群体组合
    population_combinations = [
        ('popA', 'popB'),
        ('popA', 'popC'),
        ('popA', 'popD'),
        ('popB', 'popC'),
        ('popB', 'popD'),
        ('popC', 'popD')
    ]
    
    # 创建输出目录
    output_dir = 'split_fst_files'
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理每种组合
    for pop1, pop2 in population_combinations:
        # 筛选当前组合的数据
        mask = ((df.iloc[:, 0] == pop1) & (df.iloc[:, 1] == pop2)) | \
               ((df.iloc[:, 0] == pop2) & (df.iloc[:, 1] == pop1))
        
        subset_df = df[mask].copy()
        
        if len(subset_df) == 0:
            print(f"警告: 未找到 {pop1}-{pop2} 组合的数据")
            continue
        
        # 确保群体顺序正确（字母顺序）
        swap_mask = subset_df.iloc[:, 0] != pop1
        if swap_mask.any():
            temp_col1 = subset_df.iloc[swap_mask, 0].copy()
            temp_col2 = subset_df.iloc[swap_mask, 1].copy()
            subset_df.iloc[swap_mask, 0] = temp_col2
            subset_df.iloc[swap_mask, 1] = temp_col1
        
        # 按染色体和位点排序
        subset_df = subset_df.sort_values(by=[subset_df.columns[2], subset_df.columns[3]])
        
        # 生成输出文件名
        output_filename = f"{pop1}_{pop2}_fst.txt"
        output_path = os.path.join(output_dir, output_filename)
        
        # 保存到文件
        subset_df.to_csv(output_path, sep='\t', index=False)
        print(f"已创建: {output_path} (包含 {len(subset_df)} 行数据)")

def main():
    input_file = input("请输入fst文件路径: ")
    
    if not os.path.exists(input_file):
        print(f"错误: 文件 {input_file} 不存在")
        return
    
    process_fst_file(input_file)
    print("文件处理完成！所有第六列的NA和负数值已替换为0")

if __name__ == "__main__":
    main()
