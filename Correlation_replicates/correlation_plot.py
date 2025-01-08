# 설명: 본 코드는 여러 개의 쌍으로 이루어진 데이터를 불러와서 각 쌍의 상관관계를 산점도로 그려주는 코드입니다.
# Description: Correlation analysis between replicates with scatter plots, Pearson correlation coefficients, and customized axis formatting, saving the results as a PNG file.
   
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
from PIL import Image

# 데이터를 불러오기
data_path = 'path setting'
df = pd.read_csv(f'{data_path}/input file.csv', sep=',')

# 필요한 컬럼 추출
columns = []
for col in df.columns:
    if '_' in col:
        col_splitted = col.split('_')
        if len(col_splitted) == 2:  # 조건을 만족하는 경우에만 추가
            condition, passage = col_splitted
            col = condition + '_' + passage[:-2]
            if not(col in columns):
                columns.append(col)

# 축 표시 단위 설정 함수
def axis_formatter(x, pos):
    return f'{x:.1f}'

# 그래프 생성 및 저장
num_plots_per_page = 9  # 한 페이지에 9개의 플롯
total_plots = len(columns) * 3
num_pages = (total_plots + num_plots_per_page - 1) // num_plots_per_page

# 임시 디렉토리 생성
temp_dir = os.path.join(data_path, 'temp_png')
os.makedirs(temp_dir, exist_ok=True)

for page in range(num_pages):
    fig, axs = plt.subplots(3, 3, figsize=(8, 8))  # A4 크기, 한 페이지에 9개의 플롯

    for plot_idx in range(num_plots_per_page):
        col_idx = page * num_plots_per_page + plot_idx
        if col_idx >= total_plots:
            break
        col = columns[col_idx // 3]
        rep_idx = col_idx % 3

        x = df[f'{col}R1']
        y = [df[f'{col}R2'], df[f'{col}R3'], df[f'{col}R3']][rep_idx]
        label_x = f'{col}R1'
        label_y = [f'{col}R2', f'{col}R3', f'{col}R3'][rep_idx]

        row, col_num = divmod(plot_idx, 3)
        ax = axs[row, col_num]
        ax.scatter(x, y, s=10, edgecolors='white', c='salmon')  # 점 크기 줄이기
        ax.set_xlabel(label_x, fontsize=10, fontname='Arial')  # 폰트 크기 유지
        ax.set_ylabel(label_y, fontsize=10, fontname='Arial')  # 폰트 크기 유지
        ax.xaxis.set_major_formatter(FuncFormatter(axis_formatter))
        ax.yaxis.set_major_formatter(FuncFormatter(axis_formatter))
        r, _ = pearsonr(x, y)
        ax.set_title(f'R-sq.: {r**2:.4f}', fontsize=10, fontname='Arial')  # 폰트 크기 유지

    plt.tight_layout(pad=2.0)  # 플롯 간격 늘리기

    # 각 페이지를 별도의 PNG 파일로 저장
    page_png_path = os.path.join(temp_dir, f'correlation_page_{page + 1}.png')
    fig.savefig(page_png_path, dpi=150)  # PNG에 페이지 저장, 해상도 낮추기
    plt.close(fig)

print(f"All pages saved as PNG in {temp_dir}")

# PNG 파일들을 PDF로 병합 (선택사항)
#merged_pdf_path = os.path.join(data_path, 'correlation_results_merged.pdf')
#images = [Image.open(os.path.join(temp_dir, f'corre_page_{page + 1}.png')).convert('RGB') for page in range(num_pages)]
#images[0].save(merged_pdf_path, save_all=True, append_images=images[1:])

# 임시 디렉토리 정리 (선택사항)
#import shutil
#shutil.rmtree(temp_dir)

#print(f"Merged PDF saved as {merged_pdf_path}")

