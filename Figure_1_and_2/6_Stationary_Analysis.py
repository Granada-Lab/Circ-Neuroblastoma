
"""
Plot Autocovariance. Figure 3A-C
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "cschmal.science@gmail.com"

# Import data science libraries
import pandas as pd
import numpy as np

# Import plotting libraries
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Plotting options
import seaborn as sns
palette_tab10 = sns.color_palette("tab10", 10)

# import acf functions
import statsmodels.tsa.stattools as st

# import optimization function
from scipy.optimize import curve_fit


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
reporter_subset, savename = dfB, "Bmal"
reporter_subset, savename = df, "All"

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
                    print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                elif ("MixedCluster" in cluster_assignments) and ("InfradianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "MixedCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "InfradianCluster"]["Experiment"].values
                        print(experiments_2_delete)
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
                    print(experiments_2_delete)
                    for m in experiments_2_delete:
                        #print(m)
                        df.drop(df.loc[df['Experiment']==m].index, inplace=True)
                        #pass
                elif ("MixedCluster" in cluster_assignments) and ("CircadianCluster" in cluster_assignments):
                    if sum(dftmp["ClusterAssignment"] == "CircadianCluster") > 1:
                        experiments_2_delete = dftmp.loc[dftmp["ClusterAssignment"] == "MixedCluster"]["Experiment"].values
                        print(experiments_2_delete)
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


print(df)

df_ACF = pd.read_csv("./Results/Data/ACovFAnalysis.csv", sep=",")

dfACF_B_tmp = df_ACF.loc[df_ACF["Reporter"] == "BMAL1"]
dfACF_B_tmp = dfACF_B_tmp.loc[dfACF_B_tmp["Experiment"].isin(df["Experiment"].values)]
print(Bmal1Clustering)
ClusteringBmal_tmp = [ Bmal1Clustering.query(f"Experiment=='{i}'")["ClusterAssignment"].values[0] for i in dfACF_B_tmp["Experiment"].values]
dfACF_B_tmp["ClusterAssignment"] = ClusteringBmal_tmp

dfACF_B_tmp = dfACF_B_tmp.loc[(dfACF_B_tmp["ClusterAssignment"] == "CircadianCluster") | (dfACF_B_tmp["ClusterAssignment"] == "MixedCluster") ]
# dfACF_B_tmp = dfACF_B_tmp.loc[(dfACF_B_tmp["ClusterAssignment"] == "CircadianCluster") ]

dfACF_P_tmp = df_ACF.loc[df_ACF["Reporter"] == "PER2"]
dfACF_P_tmp = dfACF_P_tmp.loc[dfACF_P_tmp["Experiment"].isin(df["Experiment"].values)]
ClusteringPer_tmp = [ Per2Clustering.query(f"Experiment=='{i}'")["ClusterAssignment"].values[0] for i in dfACF_P_tmp["Experiment"].values]
dfACF_P_tmp["ClusterAssignment"] = ClusteringPer_tmp

dfACF_P_tmp = dfACF_P_tmp.loc[(dfACF_P_tmp["ClusterAssignment"] == "CircadianCluster") | (dfACF_P_tmp["ClusterAssignment"] == "MixedCluster") ]
# dfACF_P_tmp = dfACF_P_tmp.loc[(dfACF_P_tmp["ClusterAssignment"]== "CircadianCluster")]


celllines_ = np.unique(list(dfACF_B_tmp["Cellline"].values) + list(dfACF_P_tmp["Cellline"].values) )
celllines_ = ["SKNSH", "NGP", "GIMEN", "CLBGA", "Lan5", "SKNAS", "SKNBE", "Kelly", "IMR5"]

###
# Initialize Figure
###
num_rows = 1
fig = plt.figure(figsize=(6.4*1.5, 4.8*1.25/2*num_rows))

plt.subplot(num_rows,3,2)

width = 0.35  # the width of the bars


mean_B, mean_P, std_B, std_P, names_ = [], [], [], [], []
for n, cellline_ in enumerate(celllines_):
    #print(cellline_, cellline_ in dfACF_B_tmp["Cellline"].values)
    #print(dfACF_B_tmp["Cellline"])
    if cellline_ in dfACF_B_tmp["Cellline"].values:
        values_ = dfACF_B_tmp.loc[ dfACF_B_tmp["Cellline"]==cellline_ ]["Period"].values 
        mean_B.append( np.mean( values_ ) )
        std_B.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))-width/2, values_, color=palette_tab10[0], edgecolor="k", zorder=np.inf, s=10)
    else:
        mean_B.append(np.nan)
        std_B.append(np.nan)
    if cellline_ in dfACF_P_tmp["Cellline"].values:
        values_ = dfACF_P_tmp.loc[ dfACF_P_tmp["Cellline"]==cellline_ ]["Period"].values
        mean_P.append( np.mean( values_ ) )
        std_P.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))+width/2, values_, color=palette_tab10[1], edgecolor="k", zorder=np.inf, s=10)

    else:
        mean_P.append(np.nan)
        std_P.append(np.nan)
    names_.append(cellline_)

ind = np.arange(len(names_))  # the x locations for the groups

rects1 = plt.bar(ind - width/2, mean_B, width, yerr=std_B,
                label='Bmal1')
rects2 = plt.bar(ind + width/2, mean_P, width, yerr=std_P,
                label='Per2')

plt.ylabel("Period (h)")
plt.xticks(ind, names_, rotation=90, horizontalalignment="center")
# xticklabels(names_, rotation=90, horizontalalignment="center")
# ylabel('circadian (%)')
plt.ylim(10, 35)
plt.title("E", loc="left")

print(dfACF_B_tmp)
print(dfACF_P_tmp)

print(celllines_)

###
# Amplitude or damping rate
###
plt.subplot(num_rows,3,3)

c_observable = "Amplitude"
c_observable, plot_halflife = "DampingCoeff", True


mean_B, mean_P, std_B, std_P, names_ = [], [], [], [], []
for n, cellline_ in enumerate(celllines_):
    #print(cellline_, cellline_ in dfACF_B_tmp["Cellline"].values)
    #print(dfACF_B_tmp["Cellline"])
    if cellline_ in dfACF_B_tmp["Cellline"].values:
        values_ = dfACF_B_tmp.loc[ dfACF_B_tmp["Cellline"]==cellline_ ][c_observable].values
        if c_observable == "DampingCoeff":
            if plot_halflife == True:
                values_ = np.log(2)/values_
            else:
                pass
        mean_B.append( np.mean( values_ ) )
        std_B.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))-width/2, values_, color=palette_tab10[0], edgecolor="k", zorder=np.inf, s=10)
    else:
        mean_B.append(np.nan)
        std_B.append(np.nan)
    if cellline_ in dfACF_P_tmp["Cellline"].values:
        values_ = dfACF_P_tmp.loc[ dfACF_P_tmp["Cellline"]==cellline_ ][c_observable].values
        if c_observable == "DampingCoeff":
            if plot_halflife == True:
                values_ = np.log(2)/values_
            else:
                pass
        mean_P.append( np.mean( values_ ) )
        std_P.append( np.std( values_ ) )
        plt.scatter(np.array([n]*len(values_))+width/2, values_, color=palette_tab10[1], edgecolor="k", zorder=np.inf, s=10)

    else:
        mean_P.append(np.nan)
        std_P.append(np.nan)
    names_.append(cellline_)



ind = np.arange(len(names_))  # the x locations for the groups

rects1 = plt.bar(ind - width/2, mean_B, width, yerr=std_B,
                label='Bmal1')
rects2 = plt.bar(ind + width/2, mean_P, width, yerr=std_P,
                label='Per2')

plt.ylabel(c_observable)
plt.xticks(ind, names_, rotation=90, horizontalalignment="center")
# xticklabels(names_, rotation=90, horizontalalignment="center")
plt.ylabel('Half-life of amplitude decline (h)')
# plt.ylim(10, 35)
plt.title("F", loc="left")


###
# Plot example
###

def damped_acovfunc(x, D, L, tau):
    acovf_model = D/L*np.exp(-L*x)*np.cos(2.*np.pi/tau*x)
    return acovf_model


plt.subplot(num_rows,3,1)

Bmal_DetrData = pd.read_csv("./Results/Data/Detrended_Bmal.csv", sep=";", index_col="Time")
print(Bmal_DetrData)
c_data = Bmal_DetrData["20220705_GIMEN_BMAL1_2"].values
t_detr = Bmal_DetrData.index.to_numpy()

#t = Bmal_DetrData.index
#plt.plot(t, signal)

acovf_data = st.acovf(c_data)
L_ = 0.005
#A_list.append( sqrt(out.params['D'].value/out.params['L'].value) )
D_ = max(acovf_data)**2 * L_
#print("Initial D: ", D_, max(acovf_data))
popt, pcov = curve_fit(damped_acovfunc, t_detr-t_detr[0], acovf_data, p0=[D_, L_, 24.])

plt.plot(t_detr-t_detr[0], acovf_data, label="Data", linewidth=2)
plt.plot(t_detr-t_detr[0], damped_acovfunc(t_detr-t_detr[0], popt[0], popt[1], popt[2]), linestyle="--", label="Fit", color="gray")
plt.legend(loc=0)
plt.xlabel("Time lag (d)")
#ylabel("Autocorrelation")
plt.xticks([0, 24, 48, 72, 96, 120], [0, 1, 2, 3, 4, 5])
plt.yticks([-2000, -1000, 0, 1000, 2000], ["-2000", "", "0", "", "2000"])
plt.ylabel("Autocovariance")
plt.xlim(0, 120)
plt.ylim(-2000, 2000)
plt.title("D", loc="left")


df_AC = pd.read_csv("./Results/Data/ACovFAnalysis.csv", sep=",")
dfAC_B_tmp = df_AC.loc[df_AC["Reporter"] == "BMAL1"]
dfAC_B_tmp = dfAC_B_tmp.loc[dfAC_B_tmp["Experiment"].isin(df["Experiment"].values)]
ClusteringBmal_tmp = [ Bmal1Clustering.query(f"Experiment=='{i}'")["ClusterAssignment"].values[0] for i in dfAC_B_tmp["Experiment"].values]
dfAC_B_tmp["ClusterAssignment"] = ClusteringBmal_tmp

dfAC_B_tmp = dfAC_B_tmp.loc[(dfAC_B_tmp["ClusterAssignment"] == "CircadianCluster") | (dfAC_B_tmp["ClusterAssignment"] == "MixedCluster") ]
# dfAC_B_tmp = dfAC_B_tmp.loc[(dfAC_B_tmp["ClusterAssignment"] == "CircadianCluster") ]

dfAC_P_tmp = df_AC.loc[df_AC["Reporter"] == "PER2"]
dfAC_P_tmp = dfAC_P_tmp.loc[dfAC_P_tmp["Experiment"].isin(df["Experiment"].values)]
ClusteringPer_tmp = [ Per2Clustering.query(f"Experiment=='{i}'")["ClusterAssignment"].values[0] for i in dfAC_P_tmp["Experiment"].values]
dfAC_P_tmp["ClusterAssignment"] = ClusteringPer_tmp

dfAC_P_tmp = dfAC_P_tmp.loc[(dfAC_P_tmp["ClusterAssignment"] == "CircadianCluster") | (dfAC_P_tmp["ClusterAssignment"] == "MixedCluster") ]
# dfAC_P_tmp = dfAC_P_tmp.loc[(dfAC_P_tmp["ClusterAssignment"]== "CircadianCluster")]


plt.tight_layout()

plt.savefig("./Results/Plots/Fig2.png")
plt.savefig("./Results/Plots/Fig2.svg")

plt.show()
