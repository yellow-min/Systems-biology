# 설명: KEGG pathway enrichment 분석을 수행하고, Benjamini-Hochberg 보정을 통해 p-value를 보정하는 스크립트
# 각 heatmap에는 p-value와 보정된 p-value가 포함되어 있습니다.
 
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

data_path = 'path setting'

# KEGG pathway 파일 로드
kegg_df = pd.read_csv(f'{data_path}/kegg file.csv')

# 클러스터 정보 파일 로드
cluster_df = pd.read_csv(f'{data_path}/input file.csv')

# 전체 gene 리스트
all_genes = kegg_df['Gene'].unique()
M = len(all_genes)

# 클러스터별 KEGG enrichment 분석
results = []

for cluster in cluster_df['Cluster'].unique():
    cluster_genes = cluster_df[cluster_df['Cluster'] == cluster]['Gene'].tolist()
    n = len(cluster_genes)
    
    for pathway in kegg_df['pathway'].unique():
        pathway_genes = kegg_df[kegg_df['pathway'] == pathway]['Gene'].tolist()
        k = len(pathway_genes)
        x = len(set(cluster_genes) & set(pathway_genes))
        
        # hypergeometric test 계산
        pval = hypergeom.sf(x-1, M, k, n)
        results.append({'Cluster': cluster, 'Pathway': pathway, 'Count': x, 'P-value': pval})

# 결과를 DataFrame으로 변환
results_df = pd.DataFrame(results)

# Benjamini-Hochberg 보정을 통해 corrected p-value 계산
results_df['Corrected P-value'] = multipletests(results_df['P-value'], method='fdr_bh')[1]

# 결과를 CSV 파일로 저장
results_df.to_csv(f'{data_path}/output.csv', index=False)

# 결과 확인
results_df.head()