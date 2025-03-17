
"""
Clustering of cell lines by their MRA components: Bmal1 and Per2 together. Assembly of most parts of Figure 1.
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "cschmal.science@gmail.com"

# # Set the backend before importing pyplot
# #import tkinter
import matplotlib
matplotlib.use('Agg') # Set backend here

# # Import data science libraries
# import pandas as pd
# import numpy as np

# # Import plotting libraries
import matplotlib.pyplot as plt
# #plt.switch_backend('TkAgg') #CE 29.04.2024 - keep text properties
from matplotlib.gridspec import GridSpec

# Import data science libraries
import pandas as pd
import numpy as np

#import matplotlib as mpl #CE 29.04.2024 - keep text properties
#import matplotlib as mpl
#mpl.rcParams['svg.fonttype'] = 'none'

# Import plotting libraries
#import matplotlib.pyplot as plt
#plt.switch_backend('TkAgg') #CE 29.04.2024 - keep text properties
#from matplotlib.gridspec import GridSpec

#import matplotlib as mpl #CE 29.04.2024 - keep text properties
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'


# Plotting options
import seaborn as sns
palette_tab10 = sns.color_palette("tab10", 10)

marker_dict = {'CHP212': "o", 'CLBGA':"v",  'GIMEN':"^",  'IMR5':"<",  'Kelly':">",  'Lan5':"s",  'NGP':"p",  'SKNAS':"P",  'SKNBE':"*", 'SKNSH':"X", 'SY5Y':"d", "U2OS":"D", "U2OSCRY1sKO":"x", "U2OSCRY2sKO":"+", "U2OSCRY12dKO":"4"}
cmap = plt.cm.get_cmap('tab20')
#print(cmap)

fig = plt.figure(figsize=(6.4*1.25, 4.8/2*3*1.25))


gs = GridSpec(6, 4, figure=fig)

###
# Plot exemplary detrended time series examples from good to bad
###
ax1 = fig.add_subplot(gs[0,2])
ax2 = fig.add_subplot(gs[0,3])
ax3 = fig.add_subplot(gs[1,2])
ax4 = fig.add_subplot(gs[1,3])

Bmal_DetrData = pd.read_csv("./Results/Data/Detrended_Bmal.csv", sep=";", index_col="Time")
Per2_DetrData = pd.read_csv("./Results/Data/Detrended_Per2.csv", sep=";", index_col="Time")

print(Bmal_DetrData)

#print(np.unique(Bmal_DetrData["cellline"].values))

t_detr = Bmal_DetrData.index.to_numpy()
dt = t_detr[1]-t_detr[0]
#print((range(0, int(t_detr[-1]), 24)))


ax1.plot(t_detr, Bmal_DetrData["20230815_SKNSH_BMAL1_3"].values / max(Bmal_DetrData["20230815_SKNSH_BMAL1_3"].values), label="Bmal1")
ax1.plot(t_detr, Per2_DetrData["20230808_SKNSH_PER2_3"].values / max(Per2_DetrData["20230808_SKNSH_PER2_3"].values), label="Per2")
ax1.set_xticks(list(range(0, int(t_detr[-1]), 24)), [0, 1, 2, 3, 4, 5])
ax1.xaxis.set_ticklabels([])
ax1.set_xlim(0, 120)
ax1.set_ylim(-1.15, 1.15)
#ax1.set_ylabel("normalized intensity (a.u.)")


# fig.text(0.005, 0.85, 'normalized intensity (a.u.)', va='center', rotation='vertical')
fig.text(0.5, 0.85, 'normalized intensity (a.u.)', va='center', rotation='vertical')

ax1.set_title("SKNSH", loc="left")

ax2.plot(t_detr, Bmal_DetrData["20220705_GIMEN_BMAL1_2"].values / max(Bmal_DetrData["20220705_GIMEN_BMAL1_2"].values), label="Bmal1")
ax2.plot(t_detr, Per2_DetrData["20220705_GIMEN_PER2_2"].values / max(Per2_DetrData["20220705_GIMEN_PER2_2"].values), label="Per2")
#ax2.plot(t_detr, Bmal_DetrData["20230815_SKNSH_BMAL1_3"].values)
#ax2.plot(t_detr, Per2_DetrData["20220705_GIMEN_PER2_2"].values)
#ax2.plot(t_detr, Bmal_DetrData["20230808_NGP_BMAL1_3"].values / max(Bmal_DetrData["20230808_NGP_BMAL1_3"].values))
#ax2.plot(t_detr, Per2_DetrData["20230808_NGP_PER2_2"].values / max(Per2_DetrData["20230808_NGP_PER2_2"].values))
#ax2.set_ylim(-200, 200)
ax2.set_xticks(list(range(0, int(t_detr[-1]), 24)), [0, 1, 2, 3, 4, 5])
ax2.xaxis.set_ticklabels([])
ax2.set_xlim(0, 120)
ax2.set_ylim(-1.15, 1.15)
ax2.set_title("GIMEN", loc="left")
ax2.legend(loc="upper right", fontsize=8)

ax3.plot(t_detr, Bmal_DetrData["20220705_SY5Y_BMAL1_2"].values / max(Bmal_DetrData["20220705_SY5Y_BMAL1_2"].values), label="Bmal1")
ax3.plot(t_detr, Per2_DetrData["20220705_SY5Y_PER2_2"].values / max(Per2_DetrData["20220705_SY5Y_PER2_2"].values), label="Per2")
ax3.set_xlim(0, 120)
ax3.set_xticks(list(range(0, int(t_detr[-1]), 24)), [0, 1, 2, 3, 4, 5])
ax3.set_ylim(-1.15, 1.15)
#ax3.xaxis.set_ticklabels([])
ax3.set_xlabel("time (d)")
ax3.set_title("SY5Y", loc="left")


ax4.plot(t_detr, Bmal_DetrData["20230815_CHP212_BMAL1_1"].values / max(Bmal_DetrData["20230815_CHP212_BMAL1_1"].values))
ax4.set_title("CHP212", loc="left")
ax4.set_xlim(0, 120)
ax4.set_ylim(-1.15, 1.15)
ax4.set_xticks(list(range(0, int(t_detr[-1]), 24)), [0, 1, 2, 3, 4, 5])
#ax4..set_ticklabels(list(range(0, int(t_detr[-1]), 24)))
ax4.set_xlabel("time (d)")



###
# DWT Example
###
#ax5 = fig.add_subplot(gs[0,2])
#ax5.set_title("B", loc="left")


#ax6 = fig.add_subplot(gs[1,2])
#ax6.set_title("C", loc="left")

# ax7 = fig.add_subplot(gs[0:2,2:4])
ax7 = fig.add_subplot(gs[0:2,0:2])
ax7.set_title("A", loc="left")

ax7.set_axis_off()
# # Hide X and Y axes label marks
# ax7.xaxis.set_tick_params(labelbottom=False)
# ax7.yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax7.set_xticks([])
# ax7.set_yticks([])

###
# K-means clustered PC-analysis of DWT results 
### 
#ax8 = fig.add_subplot(gs[2:4,0:2])
ax8 = fig.add_subplot(gs[4:6, 2:4])
ax8.set_title("F", loc="left")

df = pd.read_csv("./Results/Data/MRAAnalysis.csv", sep=",")

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

reporter_subset, savename = dfP, "Per"
#reporter_subset, savename = dfB, "Bmal"
reporter_subset, savename = df, "All"

sort_out_U2OS = True

###
# Filter out those cell lines that appear in both
###

Bmal1Clustering = pd.read_csv("./Results/Data/ClusteringResults_Bmal1.csv", sep=";")
Per2Clustering = pd.read_csv("./Results/Data/ClusteringResults_Per2.csv", sep=";")

#print(Bmal1Clustering)
if savename == "Bmal":
    #print(dfB)
    for i in np.unique(Bmal1Clustering["Cellline"].values):
        dftmp = Bmal1Clustering.loc[Bmal1Clustering["Cellline"] == i]
        if len(np.unique(dftmp["ClusterAssignment"].values)) > 1:
            cluster_assignments = np.unique(dftmp["ClusterAssignment"].values)
            #print(dftmp)
            if len(cluster_assignments) == 2:
                if "NoiseCluster" in cluster_assignments:
                    # delete all elements that are in the noise cluster -> remain all others
                    experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "NoiseCluster"]["Experiment"].values
                    #print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        dfB.drop(dfB.loc[dfB['Experiment']==m].index, inplace=True)
                elif ("MixedCluster" in cluster_assignments) and ("InfradianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "MixedCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "InfradianCluster"]["Experiment"].values
                        #print(experiments_2_delete)
                        for m in experiments_2_delete:
                            #print(m)
                            dfB.drop(dfB.loc[dfB['Experiment']==m].index, inplace=True)
                    else:
                        pass

elif savename == "Per":
    for j in np.unique(Per2Clustering["Cellline"].values):
        dftmp = Per2Clustering.loc[Per2Clustering["Cellline"] == j]
        if len(np.unique(dftmp["ClusterAssignment"].values)) > 1:
            cluster_assignments = np.unique(dftmp["ClusterAssignment"].values)
            #print(dftmp)
            if len(cluster_assignments) == 2:
                if "NoiseCluster" in cluster_assignments:
                    # delete all elements that are in the noise cluster -> remain all others
                    experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "NoiseCluster"]["Experiment"].values
                    #print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        dfP.drop(dfP.loc[dfP['Experiment']==m].index, inplace=True)
                        #pass
                elif ("MixedCluster" in cluster_assignments) and ("CircadianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "CircadianCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "MixedCluster"]["Experiment"].values
                        #print(experiments_2_delete)
                        for m in experiments_2_delete:
                            #print(m)
                            dfP.drop(dfP.loc[dfP['Experiment']==m].index, inplace=True)
                    else:
                        pass
            elif len(cluster_assignments) == 3:
                experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "InfradianCluster"]["Experiment"].values
                for m in experiments_2_delete:
                    #print(m)
                    dfP.drop(dfP.loc[dfP['Experiment']==m].index, inplace=True)
                #print(dftmp)

elif savename == "All":
    #print(df)
    for i in np.unique(Bmal1Clustering["Cellline"].values):
        dftmp = Bmal1Clustering.loc[Bmal1Clustering["Cellline"] == i]
        if len(np.unique(dftmp["ClusterAssignment"].values)) > 1:
            cluster_assignments = np.unique(dftmp["ClusterAssignment"].values)
            #print(dftmp)
            if len(cluster_assignments) == 2:
                if "NoiseCluster" in cluster_assignments:
                    # delete all elements that are in the noise cluster -> remain all others
                    experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "NoiseCluster"]["Experiment"].values
                    #print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                elif ("MixedCluster" in cluster_assignments) and ("InfradianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "MixedCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "InfradianCluster"]["Experiment"].values
                        #print(experiments_2_delete)
                        for m in experiments_2_delete:
                            #print(m)
                            df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                    else:
                        pass
    for j in np.unique(Per2Clustering["Cellline"].values):
        dftmp = Per2Clustering.loc[Per2Clustering["Cellline"] == j]
        if len(np.unique(dftmp["ClusterAssignment"].values)) > 1:
            cluster_assignments = np.unique(dftmp["ClusterAssignment"].values)
            #print(dftmp)
            if len(cluster_assignments) == 2:
                if "NoiseCluster" in cluster_assignments:
                    # delete all elements that are in the noise cluster -> remain all others
                    experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "NoiseCluster"]["Experiment"].values
                    #print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                        #pass
                elif ("MixedCluster" in cluster_assignments) and ("CircadianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "CircadianCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "MixedCluster"]["Experiment"].values
                        #print(experiments_2_delete)
                        for m in experiments_2_delete:
                            #print(m)
                            df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                    else:
                        pass
            elif len(cluster_assignments) == 3:
                experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "InfradianCluster"]["Experiment"].values
                for m in experiments_2_delete:
                    #print(m)
                    df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                #print(dftmp)


from sklearn.decomposition import PCA

np.random.seed(2)
#pca = PCA(n_components=4)
pca = PCA()

data4pca = reporter_subset.iloc[:, 3:7]
# data4pca = df.iloc[:, 3:7]
#data4pca = reporter_subset.iloc[:, 7:11]

print(data4pca)

data4pca.set_index(reporter_subset["Cellline"], inplace=True)

projected = pca.fit_transform(data4pca)
#print(pca.components_.T * np.sqrt(pca.explained_variance_))
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


#print(y_pred)
#print(kmeans)


###
# meshgrid plot of kmean boundaries
###

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 0.02  # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each

# Select based on actual values:
x_min, x_max = projected[:, 0].min() - 1, projected[:, 0].max() + 1
y_min, y_max = projected[:, 1].min() - 1, projected[:, 1].max() + 1

# just hard code the boundaries
x_min, x_max = -60, 80
y_min, y_max = -50, 80


xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)

c_cmap = plt.cm.Paired
c_cmap = "gray"
#plt.imshow(
    #Z,
    #interpolation="nearest",
    #extent=(xx.min(), xx.max(), yy.min(), yy.max()),
    #cmap=c_cmap,
    #aspect="auto",
    #origin="lower",
#)

#print(Z)
plt.contour(xx, yy, Z, levels=[0.99, 1.01], colors="gray")
plt.contour(xx, yy, Z, levels=[-0.01, 0.01], colors="gray")
plt.contour(xx, yy, Z, levels=[1.99, 2.01], colors="gray")
plt.contour(xx, yy, Z, levels=[2.99, 3.01], colors="gray")
#plt.contour(xx, yy, Z, color="k")

###
# Scatter plot of projected cell line behavior
###
#print(reporter_subset)
for i, j, c_c, c_cellline, rep in zip(PCA12.T[0], PCA12.T[1], y_pred, data4pca.T.keys(), reporter_subset["Experiment"]):
    #print( i, j, c_c, c_cellline, rep.split("_")[2])
    if rep.split("_")[2] == "BMAL1":
        ax8.scatter(x=i, y=j, color=palette_tab10[c_c], marker=marker_dict[c_cellline])
    elif rep.split("_")[2] == "PER2":
        ax8.scatter(x=i, y=j, color=palette_tab10[c_c], marker=marker_dict[c_cellline], facecolor="None")

for c_cellline in np.unique(data4pca.T.keys()):
    ax8.scatter([], [], marker=marker_dict[c_cellline], label=c_cellline, color="k")
ax8.legend(loc=0, fontsize=6, ncol=3)

#ax8.scatter(x=PCA12.T[0], y=PCA12.T[1], c=y_pred)
ax8.set_xlabel("principal component 1")
ax8.set_ylabel("principal component 2")
ax8.set_xlim(x_min, x_max)
ax8.set_ylim(y_min, y_max)

#print(y_pred)
#print(data4pca.index)
#print(reporter_subset["Experiment"].values)


#ax9 = fig.add_subplot(gs[2:4,2:4])
ax9 = fig.add_subplot(gs[4:6,0:2], projection='3d')

ax9.set_title("E", loc="left")

#ax9.scatter(data4pca["Circadian"], (data4pca["Infradian"]), c=y_pred)

#ax9.set_xlabel("Circadian component")
#ax9.set_ylabel("Infradian component")

c_marker = [marker_dict[m] for m in data4pca.T.keys()]
#ax9.scatter(xs=data4pca["Circadian"].values, ys=data4pca["Infradian"].values, zs=data4pca["Noise"].values, c=y_pred)

data4pca["Experiment"] = reporter_subset["Experiment"].values
data4pca["Cluster"] = y_pred
print(data4pca)

data4pca.to_csv("./Results/Data/ClusteringResults_FullDataSet.csv", sep=";")

for i, j, k, c, m, rep in zip(data4pca["Circadian"].values, data4pca["Infradian"].values, data4pca["Noise"].values, y_pred, c_marker, reporter_subset["Experiment"]):
    if rep.split("_")[2] == "BMAL1":
        ax9.scatter(xs=i, ys=j, zs=k, color=palette_tab10[c], marker=m)
    elif rep.split("_")[2] == "PER2":
        ax9.scatter(xs=i, ys=j, zs=k, color=palette_tab10[c], marker=m, facecolor="None")

ax9.plot([], [], color=palette_tab10[3], label="Infradian cluster")
ax9.plot([], [], color=palette_tab10[0], label="Mixed cluster")
ax9.plot([], [], color=palette_tab10[1], label="Circadian cluster")
ax9.plot([], [], color=palette_tab10[2], label="Noise cluster")
ax9.legend(loc="upper right", ncol=1, fontsize=7)
print("Ypred: ", y_pred)
ax9.set_xlabel("circadian (%)")
ax9.set_ylabel("infradian (%)")
ax9.set_zlabel("noise (%)")
#ax9.view_init(10, 40)
ax9.set_xlim(0, 100)
ax9.set_ylim(0, 100)
ax9.set_zlim(0, 100)


#ax10 = fig.add_subplot(gs[4:6, 0:2])
ax10 = fig.add_subplot(gs[2:4,2:4])
ax10.set_title("D", loc="left")

celllines_ = np.unique(df["Cellline"].values)

if sort_out_U2OS:
    celllines_ = ["SKNSH", "NGP", "GIMEN", "CLBGA", "Lan5", "SKNAS", "SKNBE", "SY5Y", "Kelly", "IMR5", "CHP212"]
    celllines_ticks = ["SKNSH", "NGP", "GIMEN", "CLBGA", "Lan5", "SKNAS", "SKNBE", "SY5Y", "Kelly", "IMR5", "CHP212"]
else:
    celllines_ = ["U2OS", "U2OSCRY2sKO", "SKNSH", "NGP", "GIMEN", "U2OSCRY1sKO", "CLBGA", "Lan5", "SKNAS", "U2OSCRY12dKO", "SKNBE", "SY5Y", "Kelly", "IMR5", "CHP212"]
    celllines_ticks = ["U-2 OS", "CRY2-KO", "SKNSH", "NGP", "GIMEN", "CRY1-KO", "CLBGA", "Lan5", "SKNAS", "CRY1/2-DKO", "SKNBE", "SY5Y", "Kelly", "IMR5", "CHP212"]


dfB_ = df.loc[df["Reporter"] == "BMAL1"]
dfP_ = df.loc[df["Reporter"] == "PER2"]

width = 0.35  # the width of the bars
mean_B, mean_P, std_B, std_P, names_ = [], [], [], [], []
for n, cellline_ in enumerate(celllines_):
    #print(cellline_, cellline_ in dfB_["Cellline"].values)
    #print(dfB_["Cellline"])
    if cellline_ in dfB_["Cellline"].values:
        values_ = dfB_.loc[ dfB_["Cellline"]==cellline_ ]["Circadian"].values
        mean_B.append( np.mean( values_ ) )
        std_B.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))-width/2, values_, color=palette_tab10[0], edgecolor="k", zorder=np.inf, s=10)
    else:
        mean_B.append(np.nan)
        std_B.append(np.nan)
    if cellline_ in dfP_["Cellline"].values:
        values_ = dfP_.loc[ dfP_["Cellline"]==cellline_ ]["Circadian"].values
        mean_P.append( np.mean( values_ ) )
        std_P.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))+width/2, values_, color=palette_tab10[1], edgecolor="k", zorder=np.inf, s=10)
    else:
        mean_P.append(np.nan)
        std_P.append(np.nan)
    names_.append(cellline_)
#print(dfB_)
ind = np.arange(len(names_))  # the x locations for the groups
width = 0.35  # the width of the bars

rects1 = ax10.bar(ind - width/2, mean_B, width, yerr=std_B,
                label='Bmal1')
rects2 = ax10.bar(ind + width/2, mean_P, width, yerr=std_P,
                label='Per2')
ax10.set_xticks(ind)
ax10.set_xticklabels(celllines_ticks, rotation=30, horizontalalignment="right")
ax10.set_ylabel('circadian (%)')
ax10.set_ylim(0, 100)

# hide top amd right axis spines
ax10.spines[['right', 'top']].set_visible(False)

###
# Save mean and std of
###

df_tmp = pd.DataFrame( np.array([mean_B, std_B, mean_P, std_P]).T, index=celllines_, columns=["Bmal1Mean_Circadian", "Bmal1Std_Circadian", "Per2Mean_Circadian","Per2Std_Circadian"])
#df_tmp["Cellline"] = celllines_
#df_tmp["Bmal1Mean_Circadian"] = mean_B
#df_tmp["Bmal1Std_Circadian"] = std_B
#df_tmp["Per2Mean_Circadian"] = mean_P
#df_tmp["Per2Std_Circadian"] = std_P

df_tmp.to_csv("./Results/Data/Fig1_AverageCircadianicity.csv", sep=",", index=True)
    


###
# Serves as input for Fig4
###


print(df_tmp)

###
# CWT analysis
###
# 
# Bmal_Period         = pd.read_csv("./Results/Data/CWT/Bmal_Periods.csv", sep=",", index_col=0)
# Bmal_Phase          = pd.read_csv("./Results/Data/CWT/Bmal_Phases.csv", sep=",", index_col=0)
# Bmal_Amplitude      = pd.read_csv("./Results/Data/CWT/Bmal_Amplitudes.csv", sep=",", index_col=0)
# Bmal_Power          = pd.read_csv("./Results/Data/CWT/Bmal_Power.csv", sep=",", index_col=0)
# Bmal_AverageSpec    = pd.read_csv("./Results/Data/CWT/Bmal_AverageSpec.csv", sep=",", index_col=0)
# 
# Per2_Period         = pd.read_csv("./Results/Data/CWT/Per2_Periods.csv", sep=",", index_col=0)
# Per2_Phase          = pd.read_csv("./Results/Data/CWT/Per2_Phases.csv", sep=",", index_col=0)
# Per2_Amplitude      = pd.read_csv("./Results/Data/CWT/Per2_Amplitudes.csv", sep=",", index_col=0)
# Per2_Power          = pd.read_csv("./Results/Data/CWT/Per2_Power.csv", sep=",", index_col=0)
# Per2_AverageSpec    = pd.read_csv("./Results/Data/CWT/Per2_AverageSpec.csv", sep=",", index_col=0)
# 

###
# Add placeholder axis
###
axB = fig.add_subplot(gs[2:4,0:2])
axB.set_title("C", loc="left")
axB.set_axis_off()

#print(Bmal_Period)
#print(Bmal_Phase)



#print(dfB_)

plt.tight_layout()

plt.savefig("./Results/Plots/Fig1.png")
plt.savefig("./Results/Plots/Fig1.svg")

plt.show()
