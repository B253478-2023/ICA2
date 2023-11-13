#!/usr/bin/python3

import os
import re
import sys
import urllib.parse
import urllib.request
import json
import subprocess


def ask_output_foleder(file):
    while True:
        output_folder = input(f"Please enter the output path you want to save the {file} to: ").rstrip(" /\\")
        try:
            os.makedirs(output_folder, exist_ok=True)
            print(f'The output path has been set: {output_folder}')
            print('#---------------------------')
            return output_folder
        except Exception as e:
            print(f"An error occurred while creating the folder: {e}. Please try again.")


# 获取用户想要的数据并简单分析，告知用户数据包含的序列数量和物种数，直到用户得到满意的数据
def get_data():
    print('#===========================')
    while True:
        # get input
        taxonomic_group, protein_family, number = get_user_input()
        # 搜索NCBI ID
        protein_ids = search_ncbi_ids(taxonomic_group, protein_family, number)
        # 如果检索不为空
        if protein_ids:
            #  获取FASTA数据
            fasta_data = fetch_fasta_from_ncbi(protein_ids)
        # 检索不到数据，退出程序
        else:
            print("Unable to find the data you want. Please enter the correct protein family and taxonomic group.")
            continue
        # 询问用户保存地址
        fasta_output = ask_output_foleder("fasta data")
        # 分析FASTA数据，查看包含序列数量和物种数量
        sequence_count, species_count, species_set = parse_fasta(fasta_data)
        # 告知用户序列数和物种数
        print(f"The dataset contains {sequence_count} sequences from {species_count} species.")
        # 保存fasta data来让用户检查
        fasta_file = save_fasta_to_file(fasta_data, protein_family, taxonomic_group, fasta_output)

        # 询问用户对数据是否满意，如果不满意，return none值，进入循环；如果满意，return数据，跳出循环
        if ask_continue(fasta_file):
            return fasta_file, taxonomic_group, protein_family


# 获取用户输入，想要搜索的taxonomic group和protein_family以及序列的最大数量
def get_user_input():
    # 获取用户输入的分类群
    taxonomic_group = input(
        "Please enter the taxonomic group you are interested in (e.g., Aves, Mammalia, Rodentia, Vertebrata): ")
    # 获取用户输入的蛋白质家族
    protein_family = input(
        "Please enter the protein family you want to analyze (e.g., glucose-6-phosphatase, kinases, cyclases, transporters): ")
    maxnum = input(
        "Please enter the maximum number of sequences you want to search for (Recommended: 1000; Warning: A large quantity may cause access to NCBI to fail): ")
    print(f"We will search in NCBI: {protein_family}[Title] AND {taxonomic_group}[Organism]，max number = {maxnum}")
    print('#---------------------------')
    return taxonomic_group, protein_family, maxnum


def search_ncbi_ids(taxonomic_group, protein_family, maxnum):
    # 构建Entrez esearch URL
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # 搜索条件有待商榷
    query = f"{protein_family}[Title] AND {taxonomic_group}[Organism]"
    params = dict(db="protein", term=query, retmode="json", retmax=maxnum)
    url = base_url + "?" + urllib.parse.urlencode(params)

    # 发送请求
    with urllib.request.urlopen(url) as response:
        data = response.read()
        search_results = json.loads(data)
        id_list = search_results["esearchresult"]["idlist"]
        return id_list


# 获取id对应的fasta
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


# 提取数据集中包含的序列数量和物种数量
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


# 把fasta数据保存到文件夹
def save_fasta_to_file(fasta_data, protein_family, taxonomic_group, directory):
    # 格式化文件名
    file_name = f"{protein_family}_in_{taxonomic_group}.fasta"
    # 替换可能导致文件系统问题的字符
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    # 文件路径
    fasta_path = os.path.join(directory, file_name)

    # 写入数据到文件
    with open(fasta_path, 'w') as file:
        file.write(fasta_data)
    print(f"Fasta data has been saved to: {fasta_path}. You can check it now.")

    return f"{fasta_path}"


# 询问是否要用这个数据做后续工作
def ask_continue(fasta_file):
    while True:
        # 给用户选择是否继续的选项
        choice = input(f"Do you want to continue with this dataset: {fasta_file} ? (Y/N): ").strip().upper()
        # 如果选择继续
        if choice == 'Y':
            print("Continuing with the dataset...")
            print('#===========================')
            return True
        # 如果选择换数据
        elif choice == 'N':
            print('#===========================')
            print("Please choose a different dataset.")
            return False
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")

# conservation analysis的整体模块
def conservation_analysis(protein_family, taxonomic_group, fasta_file):
    print('#===========================')
    print("We will determine and plot the level of conservation between the protein sequences！")
    con_output = ask_output_foleder('Conservation Analysis')
    #### 是否要限定序列的数量
    # 序列对齐
    aligned_file = align_sequence(protein_family, taxonomic_group, fasta_file, con_output)

    # 运行infoalign并获取输出
    infoalign_output = run_infoalign(aligned_file)
    if infoalign_output:
        # 解析输出
        sequences_info = parse_infoalign_output(infoalign_output)
        print(f"Alignment sequences have been parsed")
        # 根据同一性阈值筛选序列ID
        identity_threshold = 70
        filtered_sequence_ids = [seq_id for seq_id, Ident in sequences_info.items() if
                                 Ident >= identity_threshold]
        # 从原始FASTA文件中提取筛选出的序列
        filtered_sequences = extract_sequences(aligned_file, filtered_sequence_ids)
        print(f"Sequences with identity >= {identity_threshold}%: {filtered_sequences}")
        print(f"Alignment sequences have been filtered with identity threshold = {identity_threshold}")
        # 将筛选出的序列保存为新的FASTA文件
        selected_aligned_file = save_fasta(filtered_sequences, protein_family, taxonomic_group, con_output)

        # 可视化保守分数,需要用户输入windowsize
        plot_con(protein_family, taxonomic_group, selected_aligned_file, con_output)
    else:
        print("Error running infoalign.")

    print('#===========================')


# 序列对齐，请确保你安装了clustal omega
def align_sequence(protein_family, taxonomic_group, fasta_file, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_aligned.fasta'
    aligned_file = os.path.join(directory, file_name)

    try:
        print('Clustal Omega is working for you to align the sequences. Please be patient and wait.')
        subprocess.run(["clustalo", "-i", fasta_file, "-o", aligned_file, "--force"], check=True)
        print(f"Aligned sequences have been saved to {aligned_file}")
        return aligned_file
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")


# 调用infoalign并捕获输出
def run_infoalign(alignment_file):
    try:
        print('Infroalign is working for you to analyse the sequences.')
        cmd = ['infoalign', alignment_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print('Infoalign analyses successfully!')
        return result.stdout if result.returncode == 0 else None
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")


# 解析infoalign的输出
def parse_infoalign_output(infoalign_output):
    print("Parsing the output of infoalign.")
    # 创建一个字典来保存序列的信息
    sequences_info = {}
    # 按行分割输出文本
    lines = infoalign_output.strip().split('\n')
    # 跳过头部
    for line in lines[3:]:
        parts = line.split()
        seq_id = parts[0]
        percent_identity = float(parts[2].strip('%'))
        sequences_info[seq_id] = percent_identity
    return sequences_info


# 从FASTA文件中提取序列
def extract_sequences(fasta_file, sequence_ids):
    sequences = {}
    with open(fasta_file, 'r') as f:
        record = None
        for line in f:
            if line.startswith('>'):
                record = line.strip().split('>')[1]
                sequences[record] = ''
            elif record:
                sequences[record] += line.strip()
    # 只保留筛选出的序列
    return {seq_id: sequences[seq_id] for seq_id in sequence_ids if seq_id in sequences}


# 将序列保存为FASTA格式
def save_fasta(sequences, protein_family, taxonomic_group, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_aligned_selected.fasta'
    # 替换可能导致文件系统问题的字符
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    selected_aligned_file = os.path.join(directory, file_name)

    with open(selected_aligned_file, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f'>{seq_id}\n')
            f.write(f'{sequence}\n')
    print(f"Filtered sequences has been saved to: {selected_aligned_file}. You can check it now.")
    return selected_aligned_file


# 可视化保守水平,请确保你安装了emboss
def plot_con(protein_family, taxonomic_group, selceted_aligned_file, directory):
    #    aligned_file = f'Conservation_Analysis/{protein_family}_in_{taxonomic_group}_aligned.fasta'
    file_name = f'{protein_family}_in_{taxonomic_group}_conplot.png'
    plotcon_output = os.path.join(directory, file_name)

    try:
        print('Plotcon is working for you to plot the conservation level.')
        subprocess.run(["plotcon", "-sequence", selceted_aligned_file, "-graph", "png", "-goutfile", plotcon_output])
        print(f'Plot has been saved to {plotcon_output}')
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")


# prosite基序比对
def scan_prosite_motifs(protein_family, taxonomic_group, fasta_file):
    print('#===========================')
    print("We will scan protein sequences with motifs from the PROSITE database！")
    # 询问用户并设置输出路径
    patmatmotifs_output = ask_output_foleder("patmatmotifs")

    # 检查emboss_data的路径设置
    check_emboss_environment(fasta_file)

    # 初始化，为patmatmotifs运行做准备
    # 将序列剪成单一的片段
    sequences = cut_fasta(fasta_file)
    # 统计序列数量
    total_sequences = len(sequences)
    print(f'There are {total_sequences} sequences we are going to scan.')
    # 为打印进度条做准备
    processed_sequences = 0

    # 循环运行patmatmotifs
    for seq_id, sequence in sequences.items():
        # 调用run_patmatmotifs函数
        run_patmatmotifs(sequence, seq_id, patmatmotifs_output)
        processed_sequences += 1
        # 打印当前进度
        print_progress_bar('Patmatmotifs', processed_sequences, total_sequences)

    # 该模块运行结束
    print("All sequences have been processed.")
    print('#===========================')


# check the emboos_data path
def check_emboss_environment(fasta_file):
    emboss_data = os.environ.get('EMBOSS_DATA')
    if emboss_data is not None:
        print(f"EMBOSS_DATA: {emboss_data}")
        choice = input(f'Do you want to change the EMBOSS_DATA path? (Y/N): ').strip().upper()
        if choice == 'Y':
            emboss_data = input(f'Please enter the environment variables path: ')
            os.environ['EMBOSS_DATA'] = f'{emboss_data}'
    else:
        print("EMBOSS_DATA: Environment variables not declared")
        emboss_data = input(f'Please enter the EMBOSS_DATA path: ')
        os.environ['EMBOSS_DATA'] = f'{emboss_data}'

    # 测试用
    os.environ['EMBOSS_DATA'] = '/localdisk/home/software/EMBOSS-6.6.0/share/EMBOSS/data/'

    flag = False
    while flag == False:
        try:
            # 使用 subprocess.run 执行命令
            print(f'Checking the EMBOSS_DATA path: {emboss_data}.')
            test_output = 'test_output'
            result = subprocess.run(["patmatmotifs", "-sequence", fasta_file, "-outfile", test_output],
                                    check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)
            flag = True
        except subprocess.CalledProcessError as e:
            # 捕获错误并打印错误信息
            print(f"EMBOSS_DATA path isn't correct: {e.stderr}")
            emboss_data = input(f'Please enter the correct EMOSS_DATA path: ')
            os.environ['EMBOSS_DATA'] = f'{emboss_data}'
            flag = False
        else:
            print(f'Congratulation! It is correct!')
        finally:
            os.remove(test_output)
            print('#---------------------------')


# 切开FASTA文件，把序列单独保存
def cut_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = ''
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line.strip()[1:]
                sequence = ''
            else:
                sequence += line.strip()
        sequences[sequence_id] = sequence  # 添加最后一个序列
    return sequences


# 逐个序列运行patmatmotifs
def run_patmatmotifs(sequence, seq_id, patmatmotifs_output):
    # seq_id 包含了对序列的描述，我们只需要ID部分
    seq_id_part = seq_id.split()[0]  # 以空格分割并取第一部分

    # 创建临时文件存储单个序列
    temp_sequence_file = f'{seq_id_part}.fasta'
    with open(temp_sequence_file, 'w') as temp_file:
        temp_file.write(f'>{seq_id}\n{sequence}\n')

    # 运行patmatmotifs
    output_file = os.path.join(patmatmotifs_output, f'patmatmotifs_output_{seq_id}.txt')
    command = ['patmatmotifs', '-sequence', temp_sequence_file, '-outfile', output_file]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        print(f"\nError running patmatmotifs for sequence {seq_id}")

    os.remove(temp_sequence_file)


# 打印进度条
def print_progress_bar(progress, iteration, total, bar_length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(bar_length * iteration // total)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write(f'\r{progress}: |{bar}| {percent}% Complete'), sys.stdout.flush()
    if iteration == total:
        print()  # 打印最后的换行


def main():
    # 第一步，get data并简单分析
    fasta_file, taxonomic_group, protein_family, = get_data()
    # 仅测试用
    # taxonomic_group = f'Aves'
    # protein_family = f'glucose-6-phosphatase'
    # number = 1000

    # 进入交互界面，选择要使用的功能
    while True:
        print(f"\nWhat do you want to do next with your data in {fasta_file}?")
        print("1. Conversation analysis")
        print("2. Scan with motifs")
        print("3. Change dataset")
        print("4. Exit program")

        choice = input('Please enter the number before the option (e.g., 1, 2, 3): ')

        if choice == '1':
            # 执行 Conservation Analysis
            conservation_analysis(protein_family, taxonomic_group, fasta_file)
        elif choice == '2':
            # 执行 Scan the Prosite Motifs
            scan_prosite_motifs(protein_family, taxonomic_group, fasta_file)
        elif choice == '3':
            # 更改数据
            fasta_file, taxonomic_group, protein_family = get_data()
        elif choice == '4':
            # 退出程序
            print("Exiting the program.")
            print('#===========================')
            exit()
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")
    # 保守性分析
#    print('#===========================')
#    conservation_analysis(protein_family, taxonomic_group, fasta_file)
#    print('#===========================')
#    ask_continue()

# prosite基序比对
    scan_prosite_motifs(protein_family, taxonomic_group, fasta_file)
# 询问是否继续
    ask_continue()

if __name__ == "__main__":
    main()
