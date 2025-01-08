# 설명
#클러스터별로 구분하여 heatmap그리기
# 직접 지정할 색상 리스트
#custom_colors = ["#1A237E",  "#FFFFFF", "#F9A825"]

# Description: Generate heatmaps for each cluster using custom colors

from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# 데이터를 csv 파일로 읽어오기
data_path = 'Path setting'
data = pd.read_csv(f'{data_path}/input file.csv', index_col='Gene')

# Define the custom colormap
cmap = LinearSegmentedColormap.from_list('custom_cmap', ['#1A237E', 'white', '#F9A825'], N=256)

# Create subplots with adjusted figure size
count_cluster = list(data['Cluster'].value_counts().sort_index().values)
num_cluster = len(count_cluster)
fig, axs = plt.subplots(1, num_cluster , sharey=True, figsize=(60, 8), gridspec_kw={'width_ratios':count_cluster})  # Adjust the figsize as needed

# Plot the heatmaps
for cluster in range(num_cluster):
    cluster_data = data[data['Cluster'] == cluster][[
'condition1', 'condition2', 'condition3', 'condition4', 'condition5', 'condition6', 'condition7', 'condition8', 'condition9', 'condition10', 'condition11', 'condition12', 'condition13', 'condition14', 'condition15', 'condition16', 'condition17', 'condition18', 'condition19', 'condition20', 'condition21', 'condition22', 'condition23', 'condition24', 'condition25', 'condition26', 'condition27', 'condition28', 'condition29', 'condition30', 'condition31', 'condition32', 'condition33', 'condition34', 'condition35', 'condition36', 'condition37', 'condition38', 'condition39', 'condition40', 'condition41', 'condition42', 'condition43', 'condition44', 'condition45', 'condition46', 'condition47', 'condition48', 'condition49', 'condition50', 'condition51', 'condition52', 'condition53', 'condition54', 'condition55', 'condition56', 'condition57', 'condition58', 'condition59', 'condition60', 'condition61', 'condition62', 'condition63', 'condition64', 'condition65', 'condition66', 'condition67', 'condition68', 'condition69', 'condition70', 'condition71', 'condition72', 'condition73', 'condition74', 'condition75', 'condition76', 'condition77', 'condition78', 'condition79', 'condition80', 'condition81', 'condition82', 'condition83', 'condition84', 'condition85', 'condition86', 'condition87', 'condition88', 'condition89', 'condition90', 'condition91', 'condition92', 'condition93', 'condition94', 'condition95', 'condition96', 'condition97', 'condition98', 'condition99', 'condition100', 'condition101', 'condition102', 'condition103', 'condition104', 'condition105', 'condition106', 'condition107', 'condition108', 'condition109', 'condition110', 'condition111', 'condition112', 'condition113', 'condition114', 'condition115', 'condition116', 'condition117', 'condition118', 'condition119', 'condition120', 'condition121', 'condition122', 'condition123', 'condition124', 'condition125',
    ]]
    cluster_data_modi=cluster_data.astype(float).T
    #열개수 뽑고, shape-> 열개수 재서, <일때 continue, pass
       
    num_genes = len(cluster_data)
    #print(data.iloc[:,[0]])
    #print(data.iloc[:, 0]==0)
          
    axs[cluster].imshow(cluster_data_modi,
                        cmap=cmap,
                        vmin=-2, vmax=2, aspect='auto',
                        interpolation='nearest')

    # Set y-axis labels to gene names
    axs[cluster].set_yticks(np.arange(len(cluster_data_modi.index.tolist())))
    axs[cluster].set_yticklabels(cluster_data_modi.index.tolist())
   #print(cluster_data_modi.index.tolist())
    # input()

    # Set x-axis labels and remove tick labels
    axs[cluster].set_xlabel(f'Cluster {cluster}')

# Add colorbar
cbar_ax = fig.add_axes([1.05, 0, 0.05, 1])  # Adjust the position of the colorbar
cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-2, vmax=2)), ax=axs, orientation='vertical', fraction=0.02, pad=0.1)
cbar.set_label('Gene fitness')
    
# Save the plot as a PDF
#plt.tight_layout()
plt.savefig(f'{data_path}/output file.pdf')
# plt.show()
