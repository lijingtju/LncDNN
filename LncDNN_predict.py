import pickle
import os
import argparse
import itertools
from collections import Counter
import pandas as pd
from Bio import SeqIO
from collections import Counter
import numpy as np
import pandas as pd


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

def ENAC(fastas, window=5):
    if window < 1:
        print('Error: the sliding window should be greater than zero' + '\n\n')
        return 0
    AA = 'ACGU'
    encodings = []
    for i in fastas:
        sequence = i.strip()
        code = []
        for j in range(len(sequence)):
            if j < len(sequence) and j + window <= len(sequence):
                count = Counter(sequence[j:j + window])
                for key in count:
                    count[key] = count[key] / len(sequence[j:j + window])
                for aa in AA:
                    code.append(count[aa])
        encodings.append(code)
    num_columns = len(encodings)
    head = [f'ENAC.{i}' for i in range(num_columns)]
    pd.DataFrame(encodings).to_csv(output_file, header=head, index=False)

def ANF(fastas, output_file):
    AA = 'ACGU'
    encodings = []
    for i in fastas:
        sequence = i.strip()
        code = []
        for j in range(len(sequence)):
            code.append(sequence[0: j + 1].count(sequence[j]) / (j + 1))
        encodings.append(code)
    print(len(encodings))
    num_columns = len(encodings)
    head = [f'ANF.{i}' for i in range(num_columns)]
    pd.DataFrame(encodings).to_csv(output_file, header=head, index=False)


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

def select_SHAP_column(concanate_file, SHAP_file):
    columns_to_extract = {
        concanate_file: ['ANF.82', 'ENAC.33', 'ANF.75', 'ANF.109','ENAC.171',
                        'ENAC.164', 'ENAC.100', 'ENAC.3', 'ENAC.21','ANF.97',
                        'ANF.117','ENAC.12',
                                                     ],
    }


    # 初始化一个空的DataFrame来存放提取的数据
    combined_data = pd.DataFrame()


    base_name = os.path.basename(concanate_file)

    # 检查文件名是否在列名字典中
    if base_name in columns_to_extract:
        # 读取每个文件
        df = pd.read_csv(concanate_file)
        print(df.shape)

        # 提取指定的列
        extracted_columns = df[columns_to_extract[base_name]]
        print(extracted_columns)
        print(extracted_columns.shape)

        # 将提取的列合并到总的DataFrame中
        combined_data = pd.concat([combined_data, extracted_columns], axis=1)
        print(combined_data.shape)


    combined_data.to_csv(SHAP_file, index=False, header=True)


def main(data_path, file_name, output_path, outputfile):
    file_key_name = file_name[:-6]+""
    csv_file_name = file_key_name+ ".csv"
    fastaTOcsv(data_path + file_name, data_path + csv_file_name)
    data = pd.read_csv(data_path + csv_file_name)
    fasta = read_fasta(fasta_file)
    ENAC(fasta, path_fea+file_key_name+"_ENAC.csv")
    print('ENAC end.......')
    ANF(fasta, path_fea + file_key_name + "_ANF.csv")
    print('ANF end.......')

    ENAC_df = pd.read_csv('./features/' + file_key_name + '_ENAC.csv')
    ANF_df = pd.read_csv("./features/" + file_key_name + "_ANF.csv")

    print(ENAC_df.shape, ANF_df.shape)
    combined_features = pd.concat([ENAC_df, ANF_df], axis=1)
    concanate_file = './features/'+ path_fea+file_key_name+'ENAC_ANF.csv'
    combined_features.to_csv(concanate_file, index_label="index_label")
    SHAP_fea_file = './features/'+ path_fea+file_key_name+'ENAC_ANF_SHAP.csv'
    select_SHAP_column(concanate_file, SHAP_fea_file)
    SHAP_fea = pd.read_csv(SHAP_fea_file)

    prob, pred = get_df_prob(SHAP_fea)
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
