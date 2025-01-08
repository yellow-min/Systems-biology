# Description: Outlier detection in gene expression using Otsu thresholding and replicate treatments, with results saved as plots and CSV files.
# 설명: Otsu 임계값 설정과 반복 처리 방법을 사용하여 유전자 발현에서 이상값을 탐지하고, 결과를 그래프와 CSV 파일로 저장.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, normaltest
from skimage.filters import threshold_otsu

data_path = 'path setting'
curve=pd.DataFrame()
method = 'otsu' #실행시킬 method에 따라 바꾸기

df=pd.read_csv(f'{data_path}/input file.csv', sep=',', index_col='Gene')

def compute_threshold(data, detect_method):
    
    if detect_method=='ood':
        # OOD 
        mu = data.mean()
        sigma = data.std()
        upper_bound = mu + 2.5 * sigma
        lower_bound = mu - 2.5 * sigma
    elif detect_method=='otsu':
        # Otsu thresholding
        # https://en.wikipedia.org/wiki/Otsu's_method
        
        lower_bound = threshold_otsu(data)
        
    return lower_bound

def plot_fn(ax, condition_name, x, detect_method, replicate_treatment):
    
    data = x.values
    
    if replicate_treatment=='average':
        data = np.mean(data, axis=1) # replicate에 대해서 평균
        lower_bound = compute_threshold(data, detect_method=detect_method)
        is_outlier = (lower_bound > data)
        
        outlier_gene = df.index[is_outlier].tolist()
        num_outlier = len(outlier_gene)
        ax.set_title(f'{condition_name} (outlier : {num_outlier})')
        ax.hist(data, bins=100, color='midnightblue') #lightslategray
        ax.axvline(lower_bound, color='firebrick', linestyle='--')
        plt.savefig(f'{data_path}/outliers_{method}_{replicate_treatment}.pdf')
        
    elif replicate_treatment=='union':
        data = data.reshape(-1)
        lower_bound = compute_threshold(data, detect_method=detect_method)
        is_outlier = (lower_bound > x) # 각 replicate (열) 별로 outlier를 True / False
        is_outlier = np.any(is_outlier, axis=1) # 합집합.
        
        outlier_gene = df.index[is_outlier].tolist()
        num_outlier = len(outlier_gene)
        ax.set_title(f'{condition_name} (outlier : {num_outlier})')
        ax.hist(data, bins=100, color='midnightblue') #lightslategray
        ax.axvline(lower_bound, color='firebrick', linestyle='--')
        plt.savefig(f'{data_path}/outliers_{method}_{replicate_treatment}.pdf')
        
    elif replicate_treatment=='intersection':
        data = data.reshape(-1)
        lower_bound = compute_threshold(data, detect_method=detect_method)
        is_outlier = (lower_bound > x)
        is_outlier = np.all(is_outlier, axis=1) # 교집합.
        
        outlier_gene = df.index[is_outlier].tolist()
        num_outlier = len(outlier_gene)
        ax.set_title(f'{condition_name} (outlier : {num_outlier})')
        ax.hist(data, bins=100, color='midnightblue') #lightslategray
        ax.axvline(lower_bound, color='firebrick', linestyle='--')
        plt.savefig(f'{data_path}/outliers_{method}_{replicate_treatment}.pdf')
        
    elif replicate_treatment=='repetition':
        is_outlier_list = []
        for i in range(x.shape[1]): # 각 열에 대해서 반복
            lower_bound = compute_threshold(data[:,i], detect_method=detect_method) # 열 하나 골라서 threshold 계산
            is_outlier = (lower_bound > x.iloc[:,i])
            is_outlier_list.append(is_outlier)
            
        is_outlier = pd.concat(is_outlier_list, axis=1)
        is_outlier = np.all(is_outlier, axis=1)
        
        outlier_gene = df.index[is_outlier].tolist()
        num_outlier = len(outlier_gene)
        ax.set_title(f'{condition_name} (outlier : {num_outlier})')
        ax.hist(data[:,i], bins=100, color='midnightblue') # 마지막 replicate에 대해서 그림.
        ax.axvline(lower_bound, color='firebrick', linestyle='--')
        plt.savefig(f'{data_path}/outliers_{method}_{replicate_treatment}.pdf')
    
    return outlier_gene

condition_list = [
   'condition1', 'condition2', 'condition3', 'condition4', 'condition5', 'condition6', 'condition7', 
    ]

for replicate_treatment in ['average', 'union', 'intersection', 'repetition']:

    
    # replicate를 처리하는 방법. ['average', 'union', 'intersection', 'repetition']
    # average: replicate들을 평균값을 계산한 뒤, threshold 계산.
    # union: replicate 전체에 대해서 threshold를 한번 계산한 뒤, 합집합.
    # intersection: replicate 전체에 대해서 threshold를 한번 계산한 뒤, 교집합.
    # repetition: replicate 각각에 대해서 threshold를 계산한뒤, 교집합.    

    fig, axes = plt.subplots(4, 4, constrained_layout=True, figsize=(15,10))

    result = {}
    for i, condition_name in enumerate(condition_list):
        
        conditioned = df.loc[:, [
            f'{condition_name}1',
            f'{condition_name}2',
            f'{condition_name}3',
            f'{condition_name}4',
            f'{condition_name}5',
            f'{condition_name}6',
            f'{condition_name}7',
            ]]
        
        row = i // 4
        col = i % 4
        
        outlier_gene = plot_fn(axes[row][col], condition_name, conditioned, detect_method=method, replicate_treatment=replicate_treatment)
        result[condition_name] = outlier_gene
        
    # plt.savefig(f'{data_path}/Otsu_clustering.pdf')
    # plt.savefig(f'{data_path}/outliers_{method}_{replicate_treatment}.pdf')

    # result = {}
    # for i, condition in enumerate(df.columns):
    #      result[condition] = outlier_gene
        
    result_df = pd.DataFrame.from_dict(result, orient='index').T
    result_df.to_csv(f'{data_path}/{method}_{replicate_treatment}_outlier.csv')
