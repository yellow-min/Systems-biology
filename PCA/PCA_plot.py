# 설명: PCA 분석을 통해 유전자 데이터의 주성분을 시각화하고, 그룹별로 색상을 지정하여 결과를 PDF로 저장.
# Description: Performs PCA on gene data, visualizes principal components with group-specific colors, and saves the results as a PDF.

import pandas as pd 
import numpy as np
from sklearn.decomposition import PCA
import matplotlib as mpl
import matplotlib.pyplot as plt

# 데이터 경로 설정
data_path = 'path setting'

# 마이너스 폰트 설정
mpl.rcParams['axes.unicode_minus'] = False

# matplotlib 커스터마이징
plt.rcParams["figure.figsize"] = (20, 20)
plt.rcParams['axes.titlesize'] = 35
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['font.size'] = 30
plt.rcParams['font.family'] = 'Arial'

# 데이터 불러오기
data = pd.read_csv(f"{data_path}/input file.csv", sep=',', index_col='Gene')
data.head()

# PCA 주성분 분석
pca = PCA(random_state=1107)
X_p = pca.fit_transform(data.T)
data_pca = pd.DataFrame(X_p.T, columns=data.columns)

# 누적 설명 분산 비율
pd.Series(np.cumsum(pca.explained_variance_ratio_))

# 설명 분산 비율 퍼센트 계산
percent_variance = np.round(pca.explained_variance_ratio_ * 100, decimals=2)
columns = [f'PC{i+1}' for i in range(len(percent_variance))]

# Scree Plot 그리기
plt.bar(x=range(len(percent_variance)), height=percent_variance, tick_label=columns)
plt.ylabel('Percentage of Variance Explained')
plt.xlabel('Principal Component')
plt.title('PCA Scree Plot')
plt.show()

# 1. 특정 그룹끼리 같은 색으로 시각화 결정 (예: ab_B_3D_cecum1일 때,
# 3: '_'로 분리된 것의 앞 3개로 구분,
# 4: '_'로 분리된 것의 앞 4개로 구분)

# 컬럼 파서 함수 정의
def column_parser(col, ref_col=3):
    """
    컬럼 이름을 '_'로 분리하여 ref_col에 따라 그룹을 결정
    ref_col: 3이면 앞 3개 요소로 그룹을 결정, 4이면 앞 4개 요소로 그룹을 결정
    """
    col = col.split('_')
    if 'Pre' in col[0]:
        col = col[:2]  # 'Pre'와 그 뒤의 첫 요소까지만 남김
    elif 'F' in col[2]:
        col = col[:-1]
        col[-1] = 'F'
    else:
        if ref_col == 3:
            col = col[:-1]
        elif ref_col == 4:
            col[-1] = col[-1][:-1]
        else:
            print('ref_col is wrong!')
    col = '_'.join(col)
    return col

# 색상 사전 생성 함수 정의
def make_color_dict(data, ref_col=3):
    """
    컬럼 이름을 그룹별로 분류하고 각 그룹에 고유한 색상을 할당
    """
    color_list = []
    for col in data.columns:
        col = column_parser(col, ref_col=ref_col)
        if col not in color_list:
            color_list.append(col)

    # tab20 컬러맵 설정
    cmap = plt.cm.tab20
    color_dict = {col: cmap(float(i) / len(color_list)) for i, col in enumerate(color_list)}
    return color_dict

# 색상 사전 생성 (여기서 ref_col 값을 변경하여 앞 3개 또는 4개로 그룹화 결정)
color_dict = make_color_dict(data, 3)  # 여기서 3 또는 4로 변경 가능

# 특정 컬럼의 그룹 색상 확인
print(color_dict[column_parser(data.columns[0], 3)])  # 여기서 3 또는 4로 변경 가능

# PCA 결과 시각화
fig, axes = plt.subplots(1, 1)

for i, col in enumerate(data.columns.tolist()):
    c = column_parser(col, 3)  # 여기서 3 또는 4로 변경 가능
    axes.scatter(X_p[i, 0], X_p[i, 1], label=c, c=[color_dict[c]], alpha=0.8, s=300)

# 레전드 설정
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), fontsize=25, bbox_to_anchor=(1.05, 1), loc='upper left')

# 축 라벨 및 제목 설정
axes.set_xlabel(f'PC1: {pca.explained_variance_ratio_[0] * 100:.2f}% variance', fontsize=30)
axes.set_ylabel(f'PC2: {pca.explained_variance_ratio_[1] * 100:.2f}% variance', fontsize=30)
axes.set_title('Z75_vivo', fontsize=35)

# 레이아웃 및 저장
plt.tight_layout()
plt.savefig(f'{data_path}/output.pdf', format='pdf', bbox_inches='tight')

plt.show()
