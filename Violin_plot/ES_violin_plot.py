# Description: Violin plot visualization of values across conditions, categorized by 'Essentiality' , and saved as a PDF.
# 설명: 'Essentiality'에 따라 값을 조건별로 분류하여 바이올린 플롯으로 시각화하고 PDF로 저장.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# CSV 파일을 읽고 첫 번째 열을 인덱스로 설정
file_path = 'path setting'
df = pd.read_csv(f'{file_path}/input file.csv', index_col=0)

# 'Essentiality' 열의 데이터 타입을 문자열로 변환
df['Essentiality'] = df['Essentiality'].astype(str)

# 'Essentiality' 열에 따라 ES와 NE로 데이터를 분리
es_genes = df[df['Essentiality'] == 'ES']
ne_genes = df[df['Essentiality'] == 'NE']

# 플롯을 그릴 열 목록 (조건)
conditions = [
'condition1', 'condition2','condition3','condition4','condition5','condition6'  
    ]

# 데이터를 long format으로 변환
combined_df = df.melt(id_vars=['Essentiality'], value_vars=conditions, var_name='Condition', value_name='Value')

# 결측치 제거
combined_df = combined_df.dropna()

# 피겨 생성
plt.figure(figsize=(14, 8))

# violin plot 생성 (split 없이 hue 사용, dodge=True로 그룹 구분)
sns.violinplot(x='Condition', y='Value', hue='Essentiality', data=combined_df, inner="point", dodge=True)

# 제목 및 레이블 설정
plt.title('Violin Plot')
plt.xlabel('Condition')
plt.ylabel('Value')
plt.legend(title='Essentiality')

# 피겨를 PDF로 저장
plt.savefig(f'{file_path}/outputt.pdf')
plt.show()
