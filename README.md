# LncDNN
This repository contains machine learning models and prediction code for LncDNN. Users can utilize LncDNN to predict the localization of lncRNA (just for Nucleolus or Nucleoplasm).
![flowchart BstEP.](https://github.com/lijingtju/LncDNN/blob/main/data_process.png)

## Requirements
At the moment, a standard machine with CPUs will work.

## Installation
Currently, we suggest to run LncDNN from source.

### Manual Start
```
git clone https://github.com/lijingtju/LncDNN.git
cd /path/to/LncDNN
conda env create -f LncDNN.yml
conda activate LncDNN
mkdir features
```


### Commands to do prediction
```
python LncDNN_predict.py —fastafile ./data/test.fasta --outfile ./data/test_result.csv
```




### Note:
```
The *.fasta files should contain the sequence of lncRNA.
```
