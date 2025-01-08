# 설명: K-means 클러스터링을 수행하는 코드.
# Description: Code for performing K-means clustering

from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# 직접 지정할 색상 리스트
custom_colors = ["#1A237E",  "#FFFFFF", "#F9A825"]

# 데이터를 csv 파일로 읽어오기
data_path = 'path setting'
data = pd.read_csv(f'{data_path}/input file.csv', index_col='Gene')

# K-Means 모델 초기화
kmeans = KMeans(n_clusters=10, random_state=42, max_iter=100000) # max_iter: kmeans를 반복해서 몇번 수행할 것인가, 높을수록 성능좋아지지면 어느 숫자 이상이면 수렴함.

# 클러스터링 수행
clusters = kmeans.fit_predict(data)

# 결과 확인
result = pd.DataFrame({'Gene': data.index, 'Cluster': clusters})
#print(result)
result.to_csv(f'{data_path}/output file.csv', index=False)

# # 결과에 클러스터 정보 추가
# data['Cluster'] = clusters

# # 클러스터별로 데이터를 정렬
# sorted_data = data.sort_values(by='Cluster')

# # Create a custom colormap with three colors
# # Values below -2 are mapped to the first color, between -2 and 2 to the sequential colormap, and above 2 to the last color
# cmap_segments = [(0, '#1A237E'), (0.5, '#FFFFFF'), (1, '#F9A825')]
# custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', cmap_segments)

# # Heatmap 그리기
# plt.figure(figsize=(20, 50))
# sns.heatmap(sorted_data.drop('Cluster', axis=1), cmap=custom_cmap, cbar_kws={'shrink': 0.15}, annot=False, fmt=".2f", linewidths=0, vmax=2, vmin=-2)
# plt.title('Gene Expression Clustering Heatmap (K=15)')
# plt.savefig(f'{data_path}/K10_Clustering.png')
