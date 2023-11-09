#!/usr/bin/python3

import os
import re
import urllib.parse
import urllib.request
import json
import subprocess

#获取用户输入
def get_user_input():
    # 获取用户输入的分类群
    taxonomic_group = input(
        "Please enter the taxonomic group you are interested in (e.g., Aves, Mammalia, Rodentia, Vertebrata): ")
    # 获取用户输入的蛋白质家族
    protein_family = input(
        "Please enter the protein family you want to analyze (e.g., glucose-6-phosphatase, kinases, cyclases, transporters): ")

    print(f"You have entered Taxonomic Group: {taxonomic_group} and Protein Family: {protein_family}")

    return taxonomic_group, protein_family

#在ncbi中检索并获取id
def search_ncbi_ids(taxonomic_group, protein_family):
    # 构建Entrez esearch URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # 搜索条件有待商榷
    query = f"{protein_family}[Title] AND {taxonomic_group}[Organism]"
    params = dict(db="protein", term=query, retmode="json", retmax=1000)
    url = base_url + "?" + urllib.parse.urlencode(params)

    # 发送请求
    with urllib.request.urlopen(url) as response:
        data = response.read()
        search_results = json.loads(data)
        id_list = search_results["esearchresult"]["idlist"]
        return id_list

#获取id对应的fasta
def fetch_fasta_from_ncbi(protein_ids):
    # Construct the Entrez efetch URL for POST request
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "rettype": "fasta",
        "retmode": "text"
    }
    # Encode parameters
    data = urllib.parse.urlencode(params).encode('utf-8')
    data += b"&id=" + ",".join(protein_ids).encode('utf-8')

    # Make a POST request
    request = urllib.request.Request(base_url, data=data)
    with urllib.request.urlopen(request) as response:
        fasta_data = response.read().decode('utf-8')
        return fasta_data

#提取数据集中包含的序列数量和物种数量
def parse_fasta(fasta_data):
    # 正则表达式匹配FASTA头部，提取物种信息
    species_pattern = re.compile(r'\[([^\]]+)\]')
    species_set = set()
    sequence_count = 0

    # 分割FASTA数据到单独的序列
    sequences = fasta_data.strip().split('>')
    for seq in sequences:
        if seq:  # 忽略空字符串
            sequence_count += 1
            species_match = species_pattern.search(seq)
            if species_match:
                species = species_match.group(1)
                species_set.add(species)

    return sequence_count, len(species_set), species_set

#把fasta数据保存到文件夹
def save_fasta_to_file(fasta_data, protein_family, taxonomic_group, directory="fasta_files"):
    # 确保目录存在
    if not os.path.exists(directory):
        os.makedirs(directory)

    # 格式化文件名
    filename = f"{protein_family}_in_{taxonomic_group}.fasta"
    # 替换可能导致文件系统问题的字符
    filename = filename.replace(" ", "_").replace("/", "_").replace("\\", "_")

    # 文件路径
    file_path = os.path.join(directory, filename)

    # 写入数据到文件
    with open(file_path, 'w') as file:
        file.write(fasta_data)
    print(f"FASTA data has been saved to {file_path}")

def main():
    # get input
    taxonomic_group, protein_family = get_user_input()

    # 搜索NCBI ID
    protein_ids = search_ncbi_ids(taxonomic_group, protein_family)
    # print(protein_ids)
    if protein_ids:#如果检索不为空
        #  获取FASTA数据
        fasta_data = fetch_fasta_from_ncbi(protein_ids)
    else:#检索不到数据，退出程序
        print("Unable to find the data you want. Please enter the correct protein family and taxonomic group.")
        print("Exiting the program.")
        exit()

    # 分析FASTA数据
    sequence_count, species_count, species_set = parse_fasta(fasta_data)
    # 告知用户序列数和物种数
    print(f"The dataset contains {sequence_count} sequences from {species_count} species.")
    # 给用户选择是否继续的选项
    choice = input("Do you want to continue with this dataset? (Y/N): ").strip().upper()
    if choice == 'Y':
        print("Continuing with the dataset...")
        # 保存fasta
        save_fasta_to_file(fasta_data, protein_family, taxonomic_group)
    else:
        print("Please choose a different dataset.")
        print("Exiting the program.")
        exit()

    # 序列对齐
    fasta_file = f'fasta_files/{protein_family}_in_{taxonomic_group}.fasta'
    output_file = 'aligned_sequences.fasta'
    try:
        subprocess.run(["clustalo", "-i", fasta_file, "-o", output_file, "--force"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    main()
