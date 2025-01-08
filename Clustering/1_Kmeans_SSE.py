# 설명: K-means clustering을 이용하여 SSE를 계산하고, 그래프로 표현하는 코드.
# Description: Elbow method implementation using K-means clustering to determine the optimal number of clusters (k) and save the resulting SSE plot as a PDF.

import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# 데이터를 csv 파일로 읽어오기
data_path = 'path setting'
data = pd.read_csv(f'{data_path}/input file.csv', index_col='Gene')

# k 값에 따른 SSE를 계산하여 그래프로 표현
sse = []
for k in range(1, 26): # 확인해볼 K값
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(data)
    sse.append(kmeans.inertia_)

# SSE 그래프 표시
plt.figure(figsize=(8, 6))
plt.plot(range(1, 26), sse, marker='o') # 확인해볼 K값
plt.title('Elbow Method for Optimal k')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Sum of Squared Errors (SSE)')

# 그래프를 PDF로 저장
plt.savefig(f'{data_path}/output.pdf') # 확인해볼 K값

# 그래프 표시
#plt.show()