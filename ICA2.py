#!/usr/bin/python3

import os
import re
import sys
import urllib.parse
import urllib.request
import json
import subprocess


def ask_output_foleder(folder_name):
    while True:
        output_folder = input(f"Please enter the output path you want to save the {folder_name} to: ").rstrip(" /\\")
#       full_path = os.path.join(output_folder, folder_name)
        try:
            os.makedirs(output_folder, exist_ok=True)
        except Exception as e:
            print(f"An error occurred while creating the folder: {e}. Please try again.")
        else:
            print(f'The output path for {folder_name} has been set: {output_folder}')
            print('#---------------------------')
            return output_folder

# 打印进度条
def print_progress_bar(progress, iteration, total, bar_length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(bar_length * iteration // total)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write(f'\r{progress}: |{bar}| {percent}% Complete'), sys.stdout.flush()
    if iteration == total:
        print()  # 打印最后的换行
# 获取用户想要的数据并简单分析，告知用户数据包含的序列数量和物种数，直到用户得到满意的数据


def get_data():
    print('#===========================\n')
    while True:
        # get input
        taxonomic_group, protein_family, partial, maxnum = get_user_input()
        # 搜索NCBI ID
        fasta_data = run_edirect(taxonomic_group, protein_family, partial, maxnum)
        # 分析FASTA数据，查看包含序列数量和物种数量
        sequence_count, species_count, species_set = parse_fasta(fasta_data)
        # 告知用户序列数和物种数
        print(f"The dataset contains {sequence_count} sequences from {species_count} species.")
        fasta_output = ask_output_foleder("protein data")
        # 询问用户保存地址
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

    # 是否要保留partial序列
    while True:
        partial = input(
            "Do you want to keep the partial sequence? (Y/N): ").strip().upper()
        if partial == "Y" or partial == "N":
            break
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")

    # 获取保留序列的最大数量
    while True:
        try:
            maxnum = int(input(
                "Please enter the maximum number of sequences you want to search for (Recommended: 1000; Warning: A large quantity may cause access to NCBI to fail): "))
            if maxnum > 0:
                break
            else:
                print("Please enter a positive integer.")
        except ValueError:
            print("Invalid input. Please enter a valid positive integer.")

    if partial == "Y":
        print(f"We will search in NCBI: {protein_family}[Title] AND {taxonomic_group}[Organism]，max amount = {maxnum}...")
    elif partial == "N":
        print(f"We will search in NCBI: {protein_family}[Title] AND {taxonomic_group}[Organism] NOT partial，max amount = {maxnum}...")

    print('#---------------------------')
    return taxonomic_group, protein_family, partial, maxnum


def run_edirect(taxonomic_group, protein_family, partial, maxnum):
    if partial == "Y":
        query = f'"{protein_family}[Title] AND {taxonomic_group}[Organism]"'
    elif partial == "N":
        query = f'"{protein_family}[Title] AND {taxonomic_group}[Organism] NOT partial"'
    cmd = f'esearch -db protein -query {query} -retmax {maxnum} | efetch -format fasta'
    #efetch_cmd = 'efetch -format fasta'
    try:
        result = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        return None
    else:
        return result.stdout

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
    print("\nWe will determine and plot the level of conservation between the protein sequences！")
    con_output = ask_output_foleder('Conservation Analysis')


    # 序列对齐
    aligned_file = align_sequence(protein_family, taxonomic_group, fasta_file, 'aligned', con_output)

    # 运行infoalign并获取输出
    infoalign_output = run_infoalign(protein_family, taxonomic_group, aligned_file, con_output)

    #若infoalign有输出
    if infoalign_output:
        # 提取infoalign输出的内容
        sequences_info = parse_infoalign_output(infoalign_output)
        # 询问用户筛选条件
        filtered_sequence_ids = ask_filter(sequences_info)
        # 从原始FASTA文件中提取筛选出的序列
        filtered_sequences = extract_sequences(fasta_file, filtered_sequence_ids)
        # 将筛选出的序列保存为新的FASTA文件
        selected_fasta_file = save_selected_fasta(filtered_sequences, protein_family, taxonomic_group, con_output)
        #把筛选出的序列从小进行序列对齐
        aligned_filtered_file = align_sequence(protein_family, taxonomic_group, selected_fasta_file,'aligned_filtered', con_output)
        # 可视化保守分数,需要用户输入windowsize
        plot_con(protein_family, taxonomic_group, aligned_filtered_file, con_output)

    # 若infoalign没输出
    else:
        print("Error running infoalign.")

    # conservation analysis模块结束，打印分割线
    print('#===========================')

    return aligned_file

# 序列对齐，请确保你安装了clustal omega
def align_sequence(protein_family, taxonomic_group, fasta_file, name, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_{name}.fasta'
    aligned_file = os.path.join(directory, file_name)

    try:
        print('Clustal Omega is working for you to align the sequences. Please be patient and wait...')
        subprocess.run(["clustalo", "-i", fasta_file, "-o", aligned_file, "--force"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        print('#---------------------------')
    else:
        print(f"Aligned sequences have been saved to {aligned_file}")
        print('#---------------------------')
        return aligned_file

# 调用infoalign并捕获输出
def run_infoalign(protein_family, taxonomic_group, alignment_file, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_aligned_info'
    infoalign_file = os.path.join(directory, file_name)
    try:
        print('Infroalign is working for you to analyse the sequences...')
        subprocess.run(['infoalign', alignment_file, '-out', infoalign_file], check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        print('#---------------------------')
    else:
        print(f'Infoalign analyses result has been saved to {infoalign_file}. You can check it now.')
        print('#---------------------------')
        return infoalign_file

# 解析infoalign的输出
def parse_infoalign_output(output_file):
    print("Parsing infoalign output...")
    # 创建一个字典来保存序列的信息
    sequences_info = {}

    # 按行分割输出文本
    with open(output_file, 'r') as file:
        lines = file.readlines()

    # 跳过头部
    for line in lines[1:]:  # 第一行是标题
        parts = line.split() # id列和长度列之间有空格
        seq_id = parts[1]  # 序列ID在第二列
        parts2 = line.split('\t') # change % 后面没有空格
        percent_identity = float(parts2[-3])  # change % 在倒数第三列
        sequences_info[seq_id] = percent_identity

    print(f"Alignment sequences have been parsed.")
    print('#---------------------------')
    return sequences_info


def ask_filter(sequences_info):
    while True:
        # 根据change % 阈值筛选序列ID
        filter= input('Please enter the criteria you want to use to filter for changing values (format example: >70): ')
        if re.match(r'^>[0-9]+', filter):
            identity_threshold = float(filter.strip().split()[0][1:])
            filtered_sequence_ids = [seq_id for seq_id, percent_identity in sequences_info.items() if
                                 percent_identity >= identity_threshold]
            print('#---------------------------')
            return filtered_sequence_ids
        elif re.match(r'^<[0-9]+', filter):
            identity_threshold = float(filter.strip().split()[0][1:])
            filtered_sequence_ids = [seq_id for seq_id, percent_identity in sequences_info.items() if
                                 percent_identity <= identity_threshold]
            print('#---------------------------')
            return filtered_sequence_ids
        else:
            print("Invalid filter format. It should start with '>' or '<' followed by one or more digits!")


# 从FASTA文件中提取序列
def extract_sequences(fasta_file, filtered_sequence_ids):
    filtered_sequences = {}
    with open(fasta_file, 'r') as f:
        record = None
        for line in f:
            if line.startswith('>'):
                record = line.strip().split()[0][1:]

                if record in filtered_sequence_ids:
                    filtered_sequences[record] = ''
            elif record and record in filtered_sequence_ids:
                filtered_sequences[record] += line.strip()
    return filtered_sequences


# 将序列保存为FASTA格式
def save_selected_fasta(filtered_sequences, protein_family, taxonomic_group, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_filtered.fasta'
    # 替换可能导致文件系统问题的字符
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    filtered_file = os.path.join(directory, file_name)

    with open(filtered_file, 'w') as f:
        for seq_id, sequence in filtered_sequences.items():
            f.write(f'>{seq_id}\n')
            f.write(f'{sequence}\n')
    print(f"Filtered sequences has been saved to: {filtered_file}. You can check it now.")
    print('#---------------------------')
    return filtered_file


# 可视化保守水平,请确保你安装了emboss
def plot_con(protein_family, taxonomic_group, selceted_aligned_file, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_conplot.png'
    plotcon_output = os.path.join(directory, file_name)

    try:
        print('Plotcon is working for you to plot the conservation level.')
        subprocess.run(["plotcon", "-sequence", selceted_aligned_file, "-graph", "png", "-goutfile", plotcon_output], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    else:
        print(f'Plot has been saved to {plotcon_output}. You can check it now.')


# prosite基序比对
def scan_prosite_motifs(protein_family, taxonomic_group, fasta_file):
    print('#===========================')
    print("\nWe will scan protein sequences with motifs from the PROSITE database！")
    # 询问用户并设置输出路径
    patmatmotifs_output = ask_output_foleder("Patmatmotifs")

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
    print('#---------------------------')

    extract_motif_hits(patmatmotifs_output, protein_family, taxonomic_group)
    print('#===========================')


# check the emboos_data path
def check_emboss_environment(fasta_file):
    emboss_data = os.environ.get('EMBOSS_DATA')
    test_output = 'test_output'

    while True:
        try:
            # 使用 subprocess.run 执行命令
            print(f'Checking the EMBOSS_DATA path...')
            result = subprocess.run(["patmatmotifs", "-sequence", fasta_file, "-outfile", test_output],
                                    check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)
        except subprocess.CalledProcessError as e:
            # 捕获错误并打印错误信息
            print(f"EMBOSS_DATA path isn't correct: {e.stderr}")
            print('#---------------------------')
#            emboss_data = input(f'Please enter the correct EMOSS_DATA path: ')
#            os.environ['EMBOSS_DATA'] = f'{emboss_data}'
            if emboss_data is not None:
                print(f"Currently EMBOSS_DATA: {emboss_data}")
                emboss_data = input(f'Please enter the correct EMBOSS_DATA path: ')
                os.environ['EMBOSS_DATA'] = f'{emboss_data}'
            else:
                print("EMBOSS_DATA: Environment variables not declared")
                emboss_data = input(f'Please enter the EMBOSS_DATA path: ')
                os.environ['EMBOSS_DATA'] = f'{emboss_data}'
        else:
            print(f'Congratulation! It is correct!')
            print('#---------------------------')
            break
        finally:
            os.remove(test_output)


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

    #输出路径为：用户设置路径/patmatmotifs/patmatmotifs_output_{seq_id}.txt
    patmotmotifs_path = os.path.join(patmatmotifs_output, 'patmatmotifs')
    output_file = os.path.join(patmotmotifs_path, f'patmatmotifs_output_{seq_id}.txt')

    # 运行patmatmotifs
    command = ['patmatmotifs', '-sequence', temp_sequence_file, '-outfile', output_file]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        print(f"\nError running patmatmotifs for sequence {seq_id}")

    os.remove(temp_sequence_file)


# 提取并总结patmatmotifs输出文件中相关的motifs信息
def extract_motif_hits(output_folder, protein_family, taxonomic_group):

    # 保存路径为:用户设置路径/{protein_family}_in_{taxonomic_group}_motif_hits_summary.txt
    # 创建总结文件名
    summary_file_name = f'{protein_family}_in_{taxonomic_group}_motif_hits_summary.txt'
    summary_file = os.path.join(output_folder, summary_file_name)

    # 为打印进度条做准备
    processed_file = 0
    count = 0
    for file_name in os.listdir(output_folder):
        if file_name.startswith('patmatmotifs_output_'):
            count += 1

    # 打开总结文件
    with open(summary_file, 'w') as summary:
        # 遍历所有patmatmotifs输出文件
        for file_name in os.listdir(output_folder):
            # 确保只遍历输出文件
            if file_name.startswith('patmatmotifs_output_'):
                # 打开文件
                with open(os.path.join(output_folder, file_name), 'r') as file:
                    content = file.read()
                    # 如果hit count 不为0，即存在相关motifs
                    if "HitCount: 0" not in content:
                        # 写入文件名
                        summary.write(f"File: {file_name}\n")
                        # 从sequence行开始记录信息
                        start = content.find("# Sequence:")
                        end = content.find("#--", start)
                        summary.write(content[start:end])
                        # 分隔每个序列
                        summary.write("----------------------------------\n\n")

                processed_file += 1
                print_progress_bar('Summary', processed_file, count, bar_length=50)

    # 打印完成语句
    print(f"The relevant motif information has been saved to: {summary_file}.")


def blast_analysis(protein_family, taxonomic_group, fasta_file):
    print('#===========================')
    print(f"\nWe will use the previously obtained sequence: {fasta_file} as the blast database!！")
    print("Before conducting blast analysis, please ensure that you have obtained the sequence you want to use for blast analysis in advance! ")

    while True:
        choice = input("Do you want to continue with it? (Y/N): ").strip().upper()
            # 如果选择继续
        if choice == 'Y':
            print("Continuing!")
            print('#---------------------------')
            break
            # 如果选择换数据
        elif choice == 'N':
            print("Returning to Options Menu")
            print('#===========================')
            return
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")

    blast_output = ask_output_foleder('Blast Analysis')
    blastdb = run_mkdb(protein_family, taxonomic_group, fasta_file, blast_output)

    # 获取用户想要blast的类型
    while True:
        blast_type = input("Please enter BLAST query type (blastp, blastx): ").strip()
        if blast_type != "blastp" and blast_type != "blastx":
            print("Unable to recognize your input, please follow the prompts to enter!")
            print('#---------------------------')
        else:
            break

    # 运行 BLAST
    run_blast(protein_family, taxonomic_group, blast_type, blastdb, blast_output)
    print('#===========================')


def run_mkdb(protein_family, taxonomic_group, fasta_file, directory):
    """
    Creates a BLAST database from a FASTA file.

    Parameters:
    protein_family (str): The protein family name.
    taxonomic_group (str): The taxonomic group.
    fasta_file (str): Path to the FASTA file.
    directory (str): Directory where the BLAST database will be saved.
    """

    # 定义database名称
    file_name = f'{protein_family}_in_{taxonomic_group}_blast_db'
    # 设置路径：用户设置路径/filename
    blastdb_output = os.path.join(directory, file_name)

    try:
        print(f'Makeblastdb is working for you to make {fasta_file} as a blast database.')
        result = subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'prot', '-out', blastdb_output],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}\n{e.stderr}")
    else:
        print(f'Blast database has been saved to {blastdb_output}. You can check it now.')
        print('#---------------------------')

    return blastdb_output


def run_blast(protein_family, taxonomic_group, blast_type, db, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}_{blast_type}_result'
    blast_output = os.path.join(directory, file_name)

    while True:
        try:
            query_file = input("Please enter the path to query the file: ").strip()

            # Check if file is a FASTA file
            if not (query_file.endswith('.fasta') or query_file.endswith('.fa')) or not is_fasta(query_file):
                print("The file does not appear to be a FASTA file.")
                continue

            # Read the first sequence for type checking
            with open(query_file, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        continue  # Skip header lines
                    first_sequence = line.strip()
                    break

            # Check sequence type based on blast_type
            if blast_type == 'blastp' and not is_protein_sequence(first_sequence):
                print("The file does not contain valid protein sequences for blastp.")
                continue
            elif blast_type == 'blastn' and not is_nucleotide_sequence(first_sequence):
                print("The file does not contain valid nucleotide sequences for blastn.")
                continue

            blast_command = [blast_type, '-query', query_file, '-db', db, '-out', blast_output]
            subprocess.run(blast_command, check=True)

        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
        else:
            print(f"{blast_type} result has been saved to {blast_output}. You can check it now.")
            print('#---------------------------')
            return


def is_fasta(filename):
    with open(filename, 'r') as file:
        first_line = file.readline()
        return first_line.startswith('>')

def is_protein_sequence(sequence):
    return all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence.upper())

def is_nucleotide_sequence(sequence):
    return all(c in 'ACGTN' for c in sequence.upper())


def stats_data(protein_family, taxonomic_group, fasta_file):
    pepstats_output = ask_output_foleder('Pepstats')
    run_pepstats(protein_family, taxonomic_group, fasta_file, pepstats_output)


def run_pepstats(protein_family, taxonomic_group, fasta_file, directory):
    print('#===========================')
    print("\nWe will calculate statistics of protein properties！ e.g. Molecular weight, Number of residues, Average residue weight.")

    file_name = f'{protein_family}_in_{taxonomic_group}.pepstats'
    pepstats_output = os.path.join(directory, file_name)

    try:
        command = ["pepstats", "-sequence", fasta_file, "-outfile", pepstats_output]
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
    else:
        print(f'Pepstats result has been saved to {pepstats_output}. You can check it now.')
        print('#===========================')

def tmap(protein_family, taxonomic_group, fasta_file):
    print('#===========================')
    print("\nWe will predict and plot transmembrane segments in protein sequences.")
    tmap_output = ask_output_foleder('Tmap')
    run_tmap(protein_family, taxonomic_group, fasta_file, tmap_output)

def run_tmap(protein_family, taxonomic_group, fasta_file, directory):
    file_name = f'{protein_family}_in_{taxonomic_group}.tmap'
    tmap_output = os.path.join(directory, file_name)

    try:
        command = ['tmap', '-sequence', fasta_file, '-graph', 'png', '-outfile', tmap_output]
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
    else:
        print(f'Tmap result has been saved to {tmap_output}. You can check it now.')
        print('#===========================')


def main():
    # 第一步，获取用户输入，找到用户想要的data并简单分析
    fasta_file, taxonomic_group, protein_family, = get_data()

    # 得到data后，进入交互界面，用户选择要使用的功能，执行完功能后回到循环再次询问
    while True:
        print(f"\nWhat do you want to do next with your data in {fasta_file}?")
        print("1. Change dataset: \n\tChange the search criteria to get new data (subsequent outputs will all be based on the new dataset)")
        print("2. Conservation analysis: \n\tMultiple Sequence Alignment, Conditional Screening, Plotting Conserved Levels")
        print("3. Scan with motifs: \n\tScan protein sequence(s) of interest with motifs from the PROSITE database")
        print("4. Blast analysis: \n\tUsing the dataset as a blast database, perform blast analysis on the specified sequences")
        print("5. View statistical data: \n\tCalculate statistics of protein properties")
        print("6. Tmap: \n\tPredict and plot transmembrane segments")
        print("6. Exit program")

        choice = input('\nPlease enter the number before the option (e.g., 1, 2, 3): ')

        if choice == '1':
            # Change the search criteria to get new data (subsequent outputs will all be based on the new dataset)
            fasta_file, taxonomic_group, protein_family = get_data()
        elif choice == '2':
            # Multiple Sequence Alignment, Conditional Screening, Plotting Conserved Levels
            aligned_file = conservation_analysis(protein_family, taxonomic_group, fasta_file)
        elif choice == '3':
            # Scan protein sequence(s) of interest with motifs from the PROSITE database
            scan_prosite_motifs(protein_family, taxonomic_group, fasta_file)
        elif choice == '4':
            # Using the dataset as a blast database, perform blast analysis on the specified sequences
            blast_analysis(protein_family, taxonomic_group, fasta_file)
        elif choice == '5':
            #Calculate statistics of protein properties
            stats_data(protein_family, taxonomic_group, fasta_file)
        elif choice == '6':
            #Predict and plot transmembrane segments
            tmap(protein_family, taxonomic_group, fasta_file)
        elif choice == '7':
            # 退出程序，打印退出文本
            print("Exiting the program.")
            print('#===========================')
            exit()
        else:
            # 用户未按正确格式输入，返回循环重新输入
            print("Unable to recognize your input, please follow the prompts to enter!")


if __name__ == "__main__":
    main()
