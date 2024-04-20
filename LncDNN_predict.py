import pickle
import os
import argparse
import itertools
from collections import Counter
import pandas as pd
from Bio import SeqIO

# 定义函数将txt文件转换为csv文件
def txt_to_csv(txt_file, csv_file):
    # 读取txt文件
    df = pd.read_csv(txt_file, sep='\t', header=None)
    # 将数据写入csv文件
    df.to_csv(csv_file, index=False, header=False)

def get_df_prob(combined_features):
    prob_list = []
    models_name = "LncDNN.pkl"
    with open(model_path + models_name, 'rb') as f:
        model = pickle.load(f)
    y_proba = model.predict_proba(combined_features)[:, 1]
    y_pred = model.predict(combined_features)
    return y_proba, y_pred

def fastaTOcsv(fasta_file, csv_file):
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        header = record.description
        subcellular_localization = header.split('_')[-1]
        data.append({'Sequence': sequence, 'Subcellular_Localization': subcellular_localization})
    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)

def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def kmer(fastas, output_file):
    k = 4
    type = "RNA"
    upto = False
    normalize = True
    encoding = []
    # header = ['#', 'label']
    header = ['Sequence', 'Cytoplasm', 'Nucleus', 'Exosome', 'Chromatin']
    NA = 'ACGU'
    if type in ("DNA_features", 'RNA'):
        NA = 'ACGT'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = [name, label]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))

        for i in fastas:
            sequence = i.strip()
            kmers = kmerArray(sequence, k)
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    pd.DataFrame(encoding).to_csv(output_file, header=None, index=False)

def read_fasta(file):
    f = open(file)
    documents = f.readlines()
    string = ""
    flag = 0
    fea = []
    for document in documents:
        if document.startswith(">") and flag == 0:
            flag = 1
            continue
        elif document.startswith(">") and flag == 1:
            string = string.upper()
            fea.append(string)
            string = ""
        else:
            string += document
            string = string.strip()
            string = string.replace(" ", "")

    fea.append(string)
    f.close()
    return fea

def process_sequences(input_file, output_file, length=277):
    with open(output_file, 'w') as outfile:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)
            if len(sequence) > length:
                # 截断序列为指定长度
                processed_sequence = sequence[:length]
            elif len(sequence) < length:
                # 补充序列为指定长度
                padding = length - len(sequence)
                processed_sequence = sequence + sequence[-padding:]
            else:
                # 序列长度已经符合要求，无需处理
                processed_sequence = sequence

            # 将处理后的序列写入输出文件
            outfile.write(f">{record.id}\n")
            outfile.write(f"{processed_sequence}\n")

def main(data_path, file_name, output_path, outputfile):
    os.makedirs(current_path + "features", exist_ok=True)

    process_sequences(data_path + file_name, data_path + file_name[:-6]+"_277.fasta")
    file_key_name = file_name[:-6]+"_277"
    csv_file_name = file_key_name+ ".csv"
    fastaTOcsv(data_path + file_name, data_path + csv_file_name)
    data = pd.read_csv(data_path + csv_file_name)
    fasta = read_fasta(fasta_file)
    kmer(fasta, path_fea+file_key_name+"_kmer.csv")
    print('kmer end.......')

    fw = open(current_path+'/PseinOne2/command_RNA.sh', 'w')
    fw.write(
        'python '+current_path+'/PseinOne2/acc.py '+current_path+'/data/'+file_name+' RNA NMBAC -out '+path_fea+file_key_name+'_NMBAC.txt\n'
        'python '+current_path+'/PseinOne2/pse.py '+current_path+'/data/'+file_name+' RNA PC-PseDNC-General -out '+path_fea+file_key_name+'_PC-PseDNC-General.txt'
        )
    fw.close()
    os.system("bash "+current_path+'/features_code/PseinOne2/command_RNA.sh')
    print("bash "+current_path+'/features_code/PseinOne2/command_RNA.sh')
    print('NMBAC and PC-PseDNC-General end.......')
    txt_to_csv(path_fea + file_key_name + '_NMBAC.txt', path_fea + file_key_name + '_NMBAC.csv')
    txt_to_csv(path_fea + file_key_name + '_PC-PseDNC-General.txt', path_fea + file_key_name + '_PC-PseDNC-General.csv')

    NMBAC_df = pd.read_csv('./features/' + file_key_name + '_NMBAC.csv')
    PC_PseDNC_df = pd.read_csv('./features/' + file_key_name + '_PC-PseDNC-General.csv')
    kmer_df = pd.read_csv("./features/" + file_key_name + "_kmer.csv")

    print(kmer_df.shape, NMBAC_df.shape, PC_PseDNC_df.shape)
    combined_features = pd.concat([kmer_df, NMBAC_df, PC_PseDNC_df], axis=1)

    prob, pred = get_df_prob(combined_features)
    df_pred = pd.DataFrame(pred, columns=['predict_label'])
    df_prob = pd.DataFrame(prob, columns=['probility'])
    result_df = pd.concat([data, df_pred, df_prob], axis=1)
    result_df.to_csv(output_path + outputfile, index=False, index_label='index_label')

if __name__ == '__main__':
    current_path = os.getcwd()
    path_fea = current_path+"/features/"
    model_path = current_path+"/model/"
    parser = argparse.ArgumentParser(description="predict the localization Nucleolus or Nucleoplasm for lncRNA")
    parser.add_argument('--fastafile', action='store', dest='fastafile', required=True, \
                        help='fasta file needs including sequence')
    parser.add_argument('--outputfile', action='store', dest='outputfile', required=True, \
                        help='outputfile with save path and save name')
    args = vars(parser.parse_args())
    fasta_file = args['fastafile']
    output_file = args['outputfile']
    data_path = os.path.dirname(fasta_file) + "/"
    file_name = os.path.basename(fasta_file)
    output_path = os.path.dirname(output_file) + "/"
    outputfile = os.path.basename(output_file)
    main(data_path, file_name, output_path, outputfile)