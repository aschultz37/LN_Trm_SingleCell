'''QC for gene lists. 
   Removes genes with given prefixes (e.g., Gm, Rp).
   Uses .csv file output from Seurat.'''

import os
import numpy as np
import pandas as pd

# List of gene prefixes to remove
prefixes = ['Gm', 'Rp']

class FileExtError(Exception):
    pass 

def extract_file_tup(file_path):
    file_basename = os.path.basename(file_path)
    return os.path.splitext(file_basename)

def file_input():
    print('Enter the file path:')
    file_path = input()
    return file_path

def read_file(file_path):
    file_ext = extract_file_tup(file_path)[1]
    try:
        if file_ext == '.csv':
            return pd.read_csv(file_path)
        else:
            raise FileExtError
    except FileNotFoundError:
        raise FileNotFoundError

def output_dataframe(file_path, df):
    file_tup = extract_file_tup(file_path)
    file_ext = file_tup[1]
    outfile_path = file_tup[0] + '_trimmed' + file_ext
    if file_ext == '.csv':
        df.to_csv(outfile_path, header=True, index=False)
    else:
        raise FileExtError    

def do_trim(df):
    global prefixes
    df_t = df.T
    for col in df_t.columns:
        gene_prefix = df_t[col][0][0:2]
        print(gene_prefix)
        if gene_prefix in prefixes:
            print("pop " + str(col))
            df_t.pop(col)
            col = col-1
    return df_t

in_file = file_input()
df = read_file(in_file)
df_t_trimmed = do_trim(df)
df_trimmed = df_t_trimmed.T
output_dataframe(in_file, df_trimmed)