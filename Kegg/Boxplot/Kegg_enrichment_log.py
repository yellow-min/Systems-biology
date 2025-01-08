# Description: Boxplot visualization of  values across conditions with statistical significance annotations, using A4 landscape size.
# A4 가로 사이즈에 2개의 subplot을 그리는 코드입니다. 이 코드를 참고하여 A4 가로 사이즈에 1개의 subplot을 그리는 코드를 작성하세요.

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (11.69, 8.27)  # A4 가로 크기
plt.rcParams["pdf.fonttype"] = 42  # Ensure text is stored as text, not as shapes
plt.rcParams["font.family"] = "Arial"  # Set the font to Arial
plt.rcParams["font.size"] = 8  # 글자 크기 8 유지
plt.rcParams["lines.linewidth"] = 0.3  # 선 두께 0.3으로 설정

from collections import OrderedDict
from scipy.stats import mannwhitneyu

data_path = 'path setting'

# 데이터 파일 불러오기
df = pd.read_csv(f'{data_path}/input file.csv', index_col='Gene')
kegg_info = pd.read_csv(f'{data_path}/kegg information.csv', index_col='pathway')

category_set = list(set(kegg_info.index))
gene_set = set(df.index)

total_data = {}
# 선택한 조건 리스트
selected_conditions = ['condition 1', 'condition 2', 'condition 3', 'condition 4', 'condition 5', 'condition 6', 'condition 7']

# 데이터 처리 및 조건에 따른 데이터 분류
for cond in df.columns:
    if cond not in selected_conditions:
        continue
    cond_data = df.loc[:, f'{cond}']
    total_data[f'{cond}'] = OrderedDict()
    for category in category_set:
        category_gene = kegg_info.loc[category]['gene']
        if isinstance(category_gene, str):
            category_gene = [category_gene]
        category_gene = list(set(category_gene).intersection(gene_set))
        category_data = cond_data.loc[category_gene]
        if isinstance(category_data, pd.Series):
            category_data = category_data.to_numpy()
        total_data[f'{cond}'][f'{category}'] = category_data
    total_data[f'{cond}']['Total Gene'] = np.concatenate(list(total_data[f'{cond}'].values()))

category_set.append('Total Gene')

# 그래프 그리기
fig = plt.figure(figsize=(11.69, 8.27))
nrows = 1
ncols = 2

for idx, (cond, cond_data) in enumerate(total_data.items()):
    if idx % (nrows * ncols) == 0 and idx > 0:
        plt.tight_layout()
        plt.savefig(f'{data_path}/Sorted_Page_{idx // (nrows * ncols)}.pdf')
        plt.clf()
        fig = plt.figure(figsize=(11.69, 8.27))
    
    ax = fig.add_subplot(nrows, ncols, (idx % (nrows * ncols)) + 1)

    category = np.array(list(cond_data.keys()))
    median = [np.median(cond_data[k]) for k in category]
    median_sorted = np.argsort(median)
    category_sorted = category[median_sorted]
    cond_data = OrderedDict([(str(k), cond_data[k]) for k in category_sorted])
    
    significant_category = []
    pval_dict = {}
    for i, (k, v) in enumerate(cond_data.items()):
        if k == 'Total Gene':
            significant_category.append(k)
        else:
            stat, pval = mannwhitneyu(cond_data[k], cond_data['Total Gene'])
            if pval < 0.05:
                significant_category.append(k)
                pval_dict[k] = pval
    cond_data = OrderedDict([(k, cond_data[k]) for k in significant_category])
    
    ax.boxplot(cond_data.values(), vert=False)
    
    keys = []
    significant_idx = []
    for k in significant_category:
        if k != 'Total Gene':
            pval = pval_dict[k]
            if pval < 1e-4:
                keys.append(f'{k} (****)')
                significant_idx.append(i)
            elif pval < 1e-3:
                keys.append(f'{k} (***)')
                significant_idx.append(i)
            elif pval < 1e-2:
                keys.append(f'{k} (**)')
                significant_idx.append(i)
            elif pval < 0.05:
                keys.append(f'{k} (*)')
                significant_idx.append(i)
            else:
                keys.append(f'{k}')
        else:
            keys.append(k)
    ax.set_yticklabels(keys)
    ax.set_xscale('symlog')

    # Total Gene의 median 값 기준으로 빨간 점선 추가
    median_total_gene = np.median(cond_data['Total Gene'])
    ax.axvline(median_total_gene, color='red', linestyle='--', linewidth=0.3)

    ax.set_xlim(-10, 10)  # x축의 최대값을 10으로 설정
    ax.set_xlabel('value (log scale)')  # x축 축제목 설정
    ax.set_title(f'Condition : {cond}')
    ax.tick_params(labelsize=12)  # 축 글자 크기를 크게 설정

plt.tight_layout()
plt.savefig(f'{data_path}/out put.pdf')
plt.clf()
