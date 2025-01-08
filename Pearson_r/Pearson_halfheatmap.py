# Description: Visualization of a correlation matrix with Pearson coefficients and significance levels annotated on a heatmap, saved as a PDF and CSV file.
# 설명: 피어슨 상관 계수와 유의성 수준이 주석으로 표시된 상관 행렬을 히트맵으로 시각화하고, 결과를 PDF와 CSV 파일로 저장.

import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl

# 데이터 파일 경로 설정
data_path = 'path setting'

# Step 1: pandas를 사용하여 CSV 파일 읽기
data = pd.read_csv(f'{data_path}/input file.csv', index_col='Gene')

# Step 2: Pearson 상관 계수와 p-value 계산
correlation_matrix = data.corr()
p_values = np.zeros_like(correlation_matrix)

for i in range(len(data.columns)):
    for j in range(len(data.columns)):
        corr, p_value = pearsonr(data.iloc[:, i], data.iloc[:, j])
        p_values[i, j] = p_value

# p-value를 유의성 표기로 변환하는 함수 정의
def significance_annotation(p_value):
    if p_value <= 0.001:
        return '***'
    elif p_value <= 0.01:
        return '**'
    elif p_value <= 0.05:
        return '*'
    else:
        return 'ns'

# 유의성 수준 주석 행렬 생성
significance_annotations = np.empty_like(correlation_matrix, dtype=object)
for i in range(len(data.columns)):
    for j in range(len(data.columns)):
        significance_annotations[i, j] = significance_annotation(p_values[i, j])

# Arial 폰트 사용 설정
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42  # PDF 파일에 Type 3 폰트 사용 (편집 가능한 텍스트)
mpl.rcParams['ps.fonttype'] = 42   # PS 파일에 Type 3 폰트 사용 (편집 가능한 텍스트)

# Step 3: seaborn을 사용하여 히트맵 생성
sns.set(style="white")

# 상삼각형 마스크 생성
mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))

# matplotlib figure 설정
f, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 사이즈에 맞게 figsize 조정

# 히트맵 그리기
heatmap = sns.heatmap(correlation_matrix, mask=mask, cmap="coolwarm", vmax=1, vmin=0.6, center=0.8,
                      square=True, linewidths=0.5, annot=True, fmt='.2f', cbar_kws={"shrink": .5, "ticks": [ 0.6, 0.8, 1]}, annot_kws={"size": 10, "color": "black"})

# Colorbar 글자 크기 설정
colorbar = heatmap.collections[0].colorbar
colorbar.ax.tick_params(labelsize=10)  # Colorbar 글자 크기 조정

# 유의성 수준 주석 추가
for i in range(correlation_matrix.shape[0]):
    for j in range(correlation_matrix.shape[1]):
        if not mask[i, j] and significance_annotations[i, j] != '':
            plt.text(j + 0.5, i + 0.5, significance_annotations[i, j],
                     ha='center', va='center', color='white', fontsize=8, weight='bold')  # * size 조정

# 축 레이블 글꼴 크기 조정 및 Arial 폰트 설정
ax.set_xticklabels(ax.get_xticklabels(), size=8, rotation=45, ha="right", fontname='Arial')  # 글꼴 크기 조정
ax.set_yticklabels(ax.get_yticklabels(), size=8, rotation=0, fontname='Arial')  # 글꼴 크기 조정

# PDF 파일로 저장
plt.savefig(f'{data_path}/output file.pdf', bbox_inches='tight')

# 상관 행렬과 p-value를 CSV 파일로 저장
correlation_matrix.to_csv(f'{data_path}/output file2.csv')
p_values_df = pd.DataFrame(p_values, index=data.columns, columns=data.columns)

# 그래프 보여주기 (주석 처리됨)
# plt.show()
