# Description: Generates a heatmap visualizing -log10(corrected p-values) for KEGG enrichment analysis, with clusters and terms organized, and saves it as a PDF.
# 설명: KEGG enrichment 분석의 -log10(교정된 p-value)를 시각화한 히트맵을 생성하고, 클러스터와 용어를 정렬하여 PDF로 저장.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties

# Configure plot size and font properties
rcParams['figure.figsize'] = (15, 10)
font_prop = FontProperties(family='Arial', size=12)
plt.rcParams['pdf.fonttype'] = 42  # Ensure text is stored as text, not as paths
plt.rcParams['ps.fonttype'] = 42

data_path = 'path setting'
df = pd.read_csv(f'{data_path}/input file.csv')

# Define custom colors
colors = ['white', '#1A237E']

# Create a colormap with a gradient
cmap = LinearSegmentedColormap.from_list('CustomGradient', colors)

# Prepare the data matrix for heatmap
term_list = sorted(list(set(df['#Term'])))
term_dict = {term: i for i, term in enumerate(term_list)}

# Define the desired order of clusters
cluster_list = ['condition1', 'condition2','condition3','condition4','condition5','condition6']  # 원하는 순서로 클러스터 지정

# If you want to keep the original automatic sorting, uncomment the next line and comment out the above line
# cluster_list = sorted(list(set(df['Cluster'])))  # Automatically detect and sort cluster labels

data = np.zeros((len(term_list), len(cluster_list)))

# Fill the data matrix with -log10(corrected p-value)
for index, row in df.iterrows():
    term_idx = term_dict[row['#Term']]
    cluster_idx = cluster_list.index(row['Cluster'])  # Find the index of the cluster label
    data[term_idx, cluster_idx] = -np.log10(row['Corrected P-Value'])

# Plot the heatmap
fig, ax = plt.subplots()
im = ax.imshow(data, cmap=cmap)

# Add colorbar
cbar = ax.figure.colorbar(im, ax=ax, pad=0.02, shrink=0.5)
cbar.ax.set_ylabel('-log10(p-val)', rotation=-90, va="bottom", fontproperties=font_prop)

# Set axis labels and ticks
ax.set_xticks(np.arange(len(cluster_list)))
ax.set_xticklabels(cluster_list, fontproperties=font_prop)
ax.set_yticks(np.arange(len(term_list)))
ax.set_yticklabels(term_list, fontproperties=font_prop)

# Draw white grid lines
ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
ax.grid(which="minor", color="gray", linestyle='-', linewidth=1)
ax.tick_params(which="minor", bottom=False, left=False)

# Rotate x-axis tick labels
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

# Add text annotations in the grid
for i in range(len(term_list)):
    for j in range(len(cluster_list)):
        text = ax.text(j, i, f'{data[i, j]:.1f}', ha="center", va="center", color="w", size=8, fontproperties=font_prop)

# Set title
ax.set_title("Kegg_enrichment (corr p-val)", fontproperties=font_prop)

# Adjust margins and layout
plt.subplots_adjust(left=0.4, right=0.95, bottom=0.3, top=0.9)
fig.tight_layout(rect=[0, 0, 1, 0.96])

# Save the figure as PDF
plt.savefig(f'{data_path}/output.pdf', format='pdf')
plt.show()
