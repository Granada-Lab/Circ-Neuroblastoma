"""
Clustering of cell lines by their MRA components: Bmal1 and Per2 individually. Creation of Figure  S2.
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "cschmal.science@gmail.com"

# Import data science libraries
import pandas as pd
import numpy as np

#import matplotlib as mpl #CE 29.04.2024 - keep text properties
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'

# Import plotting libraries
import matplotlib.pyplot as plt
plt.switch_backend('Qt5Agg') #keep text properties
from matplotlib.gridspec import GridSpec

# Plotting options
import seaborn as sns
palette_tab10 = sns.color_palette("tab10", 10)

marker_dict = {'CHP212': "o", 'CLBGA':"v",  'GIMEN':"^",  'IMR5':"<",  'Kelly':">",  'Lan5':"s",  'NGP':"p",  'SKNAS':"P",  'SKNBE':"*", 'SKNSH':"X", 'SY5Y':"d"}
cmap = plt.cm.get_cmap('tab20')
print(cmap)


#Bmal_DetrData = pd.read_csv("./Results/Data/Detrended_Bmal.csv", sep=";", index_col="Time")
#Per2_DetrData = pd.read_csv("./Results/Data/Detrended_Per2.csv", sep=";", index_col="Time")

#print(np.unique(Bmal_DetrData["cellline"].values))

#t_detr = Bmal_DetrData.index.to_numpy()
#dt = t_detr[1]-t_detr[0]
#print((range(0, int(t_detr[-1]), 24)))


df = pd.read_csv("./Results/Data/MRAAnalysis.csv", sep=",")

###
# Sort out U-2 OS cell lines
###


#print(df)
#print(np.unique(list(df["Cellline"].values) ))

tmp_celline = np.unique(df["Cellline"].values)
# dict that assigns color index to cellline
Color_Dict = {}
for i, j in enumerate(tmp_celline):
    Color_Dict[j] = i
#print(Color_Dict)

# genrate data subsets based on reporter
dfB = df.loc[df["Reporter"]=="BMAL1"]
dfP = df.loc[df["Reporter"]=="PER2"]

#cmap = matplotlib.cm.get_cmap('tab20')
#print(cmap)
# reporter_subset, savename = df, "All"
# reporter_subset, savename = dfP, "Per"
reporter_subset, savename = dfB, "Bmal"

from sklearn.decomposition import PCA

np.random.seed(2)
#pca = PCA(n_components=4)
pca = PCA()

data4pca = reporter_subset.iloc[:, 3:7]
#data4pca = reporter_subset.iloc[:, 7:11]

#print(data4pca)

data4pca.set_index(reporter_subset["Cellline"], inplace=True)

projected = pca.fit_transform(data4pca)
print("Explained variance", pca.explained_variance_)
print(pca.components_.T * np.sqrt(pca.explained_variance_))
#print("Hallo, ", pca.components_.T * np.sqrt(pca.explained_variance_))
loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
loading_matrix = pd.DataFrame(loadings, columns=['PC1', 'PC2', 'PC3', 'PC4'], index=data4pca.columns)
#loading_matrix = pd.DataFrame(loadings, columns=['PC1', 'PC2', ])
#print(loading_matrix)
#for i, txt in enumerate(data4pca.T.keys()):
    ##print (i, txt)
    #plt.annotate(txt, (projected.T[0][i], projected.T[1][i]), horizontalalignment="center", verticalalignment="bottom", fontsize=8)

PCA12 = projected[:,:2]
from sklearn.cluster import KMeans
#common_params = {
    #"n_init": "auto",
    #"random_state": random_state,
#}
y_pred = KMeans(n_clusters=4).fit_predict(PCA12)

kmeans = KMeans(n_clusters=4)
kmeans.fit(PCA12)

data4pca["Experiment"] = reporter_subset["Experiment"].values
data4pca["Cluster"] = y_pred

#print(y_pred)
#print(kmeans)
#print(data4pca)

###
# Assign cluster properties
###
cluster_assignment = {}
for i in range(0, 4):
    dftmp = data4pca.loc[data4pca["Cluster"] == i]
    noise_mean = np.mean(dftmp["Noise"].values)
    ultradian_mean = np.mean(dftmp["Ultradian"].values)
    circadian_mean = np.mean(dftmp["Circadian"].values)
    infradian_mean = np.mean(dftmp["Infradian"].values)

    #print("Bmal1 cluster assignment: ", i, noise_mean, ultradian_mean, circadian_mean, infradian_mean)
    if (noise_mean > ultradian_mean) and (noise_mean > circadian_mean) and (noise_mean > infradian_mean):
        cluster_assignment[i] = "NoiseCluster"
    elif (circadian_mean > noise_mean) and (circadian_mean > ultradian_mean) and (circadian_mean > infradian_mean) and (circadian_mean > 60.):
        cluster_assignment[i] = "CircadianCluster"
    elif (infradian_mean > noise_mean) and (infradian_mean > ultradian_mean) and (infradian_mean > circadian_mean) and (infradian_mean > 60.):
        cluster_assignment[i] = "InfradianCluster"
    else:
        cluster_assignment[i] = "MixedCluster"
    #print(dftmp)

clusters_ = [cluster_assignment[i] for i in data4pca["Cluster"].values]
data4pca["ClusterAssignment"] = clusters_
#print(reporter_subset)
#print(data4pca)

data4pca.loc[data4pca.ClusterAssignment == "NoiseCluster", "Cluster"] = 0
data4pca.loc[data4pca.ClusterAssignment == "MixedCluster", "Cluster"] = 1
data4pca.loc[data4pca.ClusterAssignment == "CircadianCluster", "Cluster"] = 2
data4pca.loc[data4pca.ClusterAssignment == "InfradianCluster", "Cluster"] = 3

data4pca.to_csv("./Results/Data/ClusteringResults_Bmal1.csv", sep=";")

fig = plt.figure(figsize=(6.4/2*3, 4.8/2*2.25) )
###
# Scatter plot of projected cell line behavior
###
plt.subplot(231)
print(cmap)
for i, j, c_c, c_cellline in zip(PCA12.T[0], PCA12.T[1], data4pca["Cluster"].values, data4pca.T.keys()):
    plt.scatter(x=i, y=j, color=palette_tab10[c_c], marker=marker_dict[c_cellline])

for c_cellline in np.unique(data4pca.T.keys()):
    plt.scatter([], [], marker=marker_dict[c_cellline], label=c_cellline, color="k")
plt.legend(loc=0, fontsize=5, ncol=2)

#ax8.scatter(x=PCA12.T[0], y=PCA12.T[1], c=y_pred)
plt.xlabel("principal component 1")
plt.ylabel("principal component 2")
#print(y_pred)
#print(data4pca.index)
#print(reporter_subset["Experiment"].values)
plt.title("A", loc="left")

plt.subplot(232)
#ax9.set_title("F", loc="left")

pc = plt.pcolor(loading_matrix, edgecolors="k", linewidth=4, cmap="PiYG")

for ii,i in enumerate(loading_matrix.index):
    for jj,j in enumerate(loading_matrix.keys()):
        #print(ii, i, jj, j)
        plt.text(ii+0.5,jj+0.5, round(loading_matrix.to_numpy()[jj][ii]), verticalalignment="center", horizontalalignment="center")
plt.xticks(np.arange(0, loading_matrix.shape[0], 1)+0.5, loading_matrix.columns)
plt.yticks(np.arange(0, loading_matrix.shape[1], 1)+0.5, loading_matrix.index)
#plt.scatter(data4pca["Circadian"], (data4pca["Infradian"]), c=y_pred)
print(loading_matrix)
#plt.xlabel("Circadian component")
#plt.ylabel("Infradian component")

plt.title("B", loc="left")


#plt.subplot(233, projection="3d")
ax = fig.add_subplot(233,projection='3d')
c_marker = [marker_dict[m] for m in data4pca.T.keys()]
#ax.scatter(xs=data4pca["Circadian"].values, ys=data4pca["Infradian"].values, zs=data4pca["Noise"].values, c=y_pred)
for i, j, k, c, m in zip(data4pca["Circadian"].values, data4pca["Infradian"].values, data4pca["Noise"].values, data4pca["Cluster"].values, c_marker):
    ax.scatter(xs=i, ys=j, zs=k, color=palette_tab10[c], marker=m)
ax.set_xlabel("circadian (%)")
ax.set_ylabel("infradian (%)")
ax.set_zlabel("noise (%)")
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_zlim(0, 100)

#print(data4pca)
#print(c_marker)
ax.set_title("C", loc="left")

reporter_subset, savename = dfP, "Per"
#reporter_subset, savename = dfB, "Bmal"
#reporter_subset, savename = df, "All"

from sklearn.decomposition import PCA

np.random.seed(2)
#pca = PCA(n_components=4)
pca = PCA()

data4pca = reporter_subset.iloc[:, 3:7]
#data4pca = reporter_subset.iloc[:, 7:11]

#print(data4pca)

data4pca.set_index(reporter_subset["Cellline"], inplace=True)

projected = pca.fit_transform(data4pca)
#print(pca.components_.T * np.sqrt(pca.explained_variance_))
#print("Hallo, ", pca.components_.T * np.sqrt(pca.explained_variance_))
loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
loading_matrix = pd.DataFrame(loadings, columns=['PC1', 'PC2', 'PC3', 'PC4'], index=data4pca.columns)
#loading_matrix = pd.DataFrame(loadings, columns=['PC1', 'PC2', ])
print(loading_matrix)
#for i, txt in enumerate(data4pca.T.keys()):
    ##print (i, txt)
    #plt.annotate(txt, (projected.T[0][i], projected.T[1][i]), horizontalalignment="center", verticalalignment="bottom", fontsize=8)

PCA12 = projected[:,:2]
from sklearn.cluster import KMeans
#common_params = {
    #"n_init": "auto",
    #"random_state": random_state,
#}
y_pred = KMeans(n_clusters=4).fit_predict(PCA12)

kmeans = KMeans(n_clusters=4)
kmeans.fit(PCA12)

print(y_pred)
print(kmeans)


data4pca["Experiment"] = reporter_subset["Experiment"].values
data4pca["Cluster"] = y_pred

print(y_pred)
print(kmeans)
print(data4pca)

###
# Assign cluster properties
###
cluster_assignment = {}
for i in range(0, 4):
    dftmp = data4pca.loc[data4pca["Cluster"] == i]
    noise_mean = np.mean(dftmp["Noise"].values)
    ultradian_mean = np.mean(dftmp["Ultradian"].values)
    circadian_mean = np.mean(dftmp["Circadian"].values)
    infradian_mean = np.mean(dftmp["Infradian"].values)

    print("Per2 cluster assignment: ", i, noise_mean, ultradian_mean, circadian_mean, infradian_mean)
    if (noise_mean > ultradian_mean) and (noise_mean > circadian_mean) and (noise_mean > infradian_mean):
        cluster_assignment[i] = "NoiseCluster"
    elif (circadian_mean > noise_mean) and (circadian_mean > ultradian_mean) and (circadian_mean > infradian_mean) and (circadian_mean > 60.):
        cluster_assignment[i] = "CircadianCluster"
    elif (infradian_mean > noise_mean) and (infradian_mean > ultradian_mean) and (infradian_mean > circadian_mean) and (infradian_mean > 60.):
        cluster_assignment[i] = "InfradianCluster"
    else:
        cluster_assignment[i] = "MixedCluster"
    #print(dftmp)

clusters_ = [cluster_assignment[i] for i in data4pca["Cluster"].values]
data4pca["ClusterAssignment"] = clusters_
#print(reporter_subset)

data4pca.loc[data4pca.ClusterAssignment == "NoiseCluster", "Cluster"] = 0
data4pca.loc[data4pca.ClusterAssignment == "MixedCluster", "Cluster"] = 1
data4pca.loc[data4pca.ClusterAssignment == "CircadianCluster", "Cluster"] = 2
data4pca.loc[data4pca.ClusterAssignment == "InfradianCluster", "Cluster"] = 3

print(data4pca)

data4pca.to_csv("./Results/Data/ClusteringResults_Per2.csv", sep=";")



###
# Scatter plot of projected cell line behavior
###
plt.subplot(234)
print(cmap)
for i, j, c_c, c_cellline in zip(PCA12.T[0], PCA12.T[1], data4pca["Cluster"].values, data4pca.T.keys()):
    plt.scatter(x=i, y=j, color=palette_tab10[c_c], marker=marker_dict[c_cellline])

for c_cellline in np.unique(data4pca.T.keys()):
    plt.scatter([], [], marker=marker_dict[c_cellline], label=c_cellline, color="k")
#plt.legend(loc=0, fontsize=8)

#ax8.scatter(x=PCA12.T[0], y=PCA12.T[1], c=y_pred)
plt.xlabel("principal component 1")
plt.ylabel("principal component 2")
#print(y_pred)
#print(data4pca.index)
#print(reporter_subset["Experiment"].values)
plt.title("D", loc="left")

plt.subplot(235)
#ax9.set_title("F", loc="left")

pc = plt.pcolor(loading_matrix, edgecolors="k", linewidth=4, cmap="PiYG")

for ii,i in enumerate(loading_matrix.index):
    for jj,j in enumerate(loading_matrix.keys()):
        #print(ii, i, jj, j)
        plt.text(ii+0.5,jj+0.5, round(loading_matrix.to_numpy()[jj][ii]), verticalalignment="center", horizontalalignment="center")
plt.xticks(np.arange(0, loading_matrix.shape[0], 1)+0.5, loading_matrix.columns)
plt.yticks(np.arange(0, loading_matrix.shape[1], 1)+0.5, loading_matrix.index)
#plt.scatter(data4pca["Circadian"], (data4pca["Infradian"]), c=y_pred)
print(loading_matrix)
#plt.xlabel("Circadian component")
#plt.ylabel("Infradian component")

plt.title("E", loc="left")


#plt.subplot(233, projection="3d")
ax = fig.add_subplot(236,projection='3d')
c_marker = [marker_dict[m] for m in data4pca.T.keys()]
#ax.scatter(xs=data4pca["Circadian"].values, ys=data4pca["Infradian"].values, zs=data4pca["Noise"].values, c=y_pred)
for i, j, k, c, m in zip(data4pca["Circadian"].values, data4pca["Infradian"].values, data4pca["Noise"].values, data4pca["Cluster"].values, c_marker):
    ax.scatter(xs=i, ys=j, zs=k, color=palette_tab10[c], marker=m)
ax.set_xlabel("circadian (%)")
ax.set_ylabel("infradian (%)")
ax.set_zlabel("noise (%)")
#ax.view_init(10, 40)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_zlim(0, 100)

#print(data4pca)
#print(c_marker)
ax.set_title("F", loc="left")

## The fix
#for spine in ax.spines.values():
    #spine.set_visible(False)

plt.tight_layout()

plt.savefig("./Results/Plots/SupFigure_SortOutCellLines.png")
plt.savefig("./Results/Plots/SupFigure_SortOutCellLines.svg")

plt.show()
