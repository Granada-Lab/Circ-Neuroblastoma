#!/usr/bin/env python

"""
Analyze Bmal1 and Per2 neuroblastom cell-line time series for their rhythmic properties.
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "cschmal.science@gmail.com"


# Plotting libraries
from pylab import*

# General libraries
import numpy as np

# Data wrangling libraries
import pandas as pd

# Detrending functions
from pyboat import sinc_smooth

# Time series analysis functions
from astropy.timeseries import LombScargle

# wavelet libraries
import pywt
# we make use of the modwt.py library from https://github.com/pistonly/modwtpy
from modwt import*

#import matplotlib as mpl #CE 29.04.2024 - keep text properties
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'

# Import plotting libraries
import matplotlib.pyplot as plt
plt.switch_backend('Qt5Agg') #keep text properties

# define reporter names
BmalRep = "Bmal1-Luc"
Per2Rep = "Per2::LUC"


# set up detrending parameters
dt          = 10./60     # in [h]
T_cut_off   = 48  # cut-off period for sinc-filter
skip        = int(5 // dt)


Bmal_DetrData = pd.read_csv("./Results/Data/Detrended_Bmal.csv", sep=";", index_col="Time")
Per2_DetrData = pd.read_csv("./Results/Data/Detrended_Per2.csv", sep=";", index_col="Time")

t_detr = Bmal_DetrData.index.to_numpy()


# Bmal1
Bmal_Columns = Bmal_DetrData.columns.to_list()
Bmal_Celline = np.array([i.split(sep="_")[1] for i in Bmal_Columns])

# Per2
Per2_Columns = Per2_DetrData.columns.to_list()
Per2_Celline = np.array([i.split(sep="_")[1] for i in Per2_Columns])

Celllines = np.unique( np.unique(Bmal_Celline).tolist() + np.unique(Per2_Celline).tolist() )


###
# Multiresolution analysis (MRA)
###

def PlotMRA_CoarseGrained(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=12.):
    f, axarr = plt.subplots(5, sharex=True, figsize=(6/1.7, 12/2.5))
    c_max = ceil(max(signal.flatten()))

    i_edge = round(t_edge/dt)

    OverallEnergy = sum([sum(n[i_edge:-i_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-i_edge])
    D1Energy = sum(MRA[0][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy


    D12  = MRA[0] + MRA[1]              # 1-2, 2-4
    D34  = MRA[2] + MRA[3]              # 4-8, 8-16
    D5   = MRA[4]                       # 16-32
    D678 = MRA[5] + MRA[6] + MRA[7]     # 32 - 64, 64 - 128, Smooth

    c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy+D8Energy]
    for i, j in enumerate([D12, D34, D5, D678]):
        axarr[i].plot(t, j, color="gray", linestyle="--")

        c_maxtmp = max(abs(j[i_edge:-i_edge]))
        axarr[i].plot(t[i_edge:-i_edge], j[i_edge:-i_edge], label = str(round(c_energy[i]*100,2)))
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylim(-c_maxtmp, c_maxtmp)
        axarr[i].set_ylabel(c_label[i])
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)


    axarr[4].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[4].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[4].set_ylim(-c_max, c_max)
    axarr[4].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[4].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[4].set_xlabel("Time in h")
    f.tight_layout()


def Calculate_CoarseGrained_MRA_coefficients(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=12.):
    c_max = ceil(max(signal.flatten()))

    i_edge = round(t_edge/dt)

    OverallEnergy = sum([sum(n[i_edge:-i_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-i_edge])
    D1Energy = sum(MRA[0][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy

    ## computes the signal summed over multiple details
    #D12  = MRA[0] + MRA[1]              # 1-2, 2-4
    #D34  = MRA[2] + MRA[3]              # 4-8, 8-16
    #D5   = MRA[4]                       # 16-32
    #D678 = MRA[5] + MRA[6] + MRA[7]     # 32 - 64, 64 - 128, Smooth

    #c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy]
    c_energy = np.array(c_energy)*100
    return c_energy

def PlotMRA_CoarseGrained_split(t, signal, MRA, dt=30./60, NumLevels=7, l_edge=12., r_edge=12.):
    f, axarr = plt.subplots(5, sharex=True, figsize=(6/1.7, 12/2.5))
    c_max = ceil(max(signal.flatten()))

    il_edge = round(l_edge/dt)
    ir_edge = round(r_edge/dt)

    OverallEnergy = sum([sum(n[il_edge:-ir_edge]**2) for n in MRA[:-1]])/len(t[il_edge:-ir_edge])
    D1Energy = sum(MRA[0][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy


    D12  = MRA[0] + MRA[1]              # 1-2, 2-4
    D34  = MRA[2] + MRA[3]              # 4-8, 8-16
    D5   = MRA[4]                       # 16-32
    D678 = MRA[5] + MRA[6] + MRA[7]     # 32 - 64, 64 - 128, Smooth

    c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy+D8Energy]
    for i, j in enumerate([D12, D34, D5, D678]):
        axarr[i].plot(t, j, color="gray", linestyle="--")

        c_maxtmp = max(abs(j[il_edge:-ir_edge]))
        axarr[i].plot(t[il_edge:-ir_edge], j[il_edge:-ir_edge], label = str(round(c_energy[i]*100,2)) + "%")
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylim(-c_maxtmp, c_maxtmp)
        axarr[i].set_ylabel(c_label[i])
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)


    axarr[4].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[4].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[4].set_ylim(-c_max, c_max)
    axarr[4].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[4].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[4].set_xlabel("Time in h")
    f.tight_layout()


def Calculate_CoarseGrained_MRA_coefficients_split(t, signal, MRA, dt=30./60, NumLevels=7, l_edge=12., r_edge=12.):
    c_max = ceil(max(signal.flatten()))

    il_edge = round(l_edge/dt)
    ir_edge = round(r_edge/dt)

    OverallEnergy = sum([sum(n[il_edge:-ir_edge]**2) for n in MRA[:-1]])/len(t[il_edge:-ir_edge])
    D1Energy = sum(MRA[0][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][il_edge:-ir_edge]**2) / len(t[il_edge:-ir_edge]) / OverallEnergy

    ## computes the signal summed over multiple details
    #D12  = MRA[0] + MRA[1]              # 1-2, 2-4
    #D34  = MRA[2] + MRA[3]              # 4-8, 8-16
    #D5   = MRA[4]                       # 16-32
    #D678 = MRA[5] + MRA[6] + MRA[7]     # 32 - 64, 64 - 128, Smooth

    #c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy]
    c_energy = np.array(c_energy)*100
    return c_energy

### 
# Make example calculation
###

# Downsample data to obtain a proper circadian frequency band
t_detr = t_detr[::3]
signal = Bmal_DetrData["20220705_GIMEN_BMAL1_2"].to_numpy()[::3]
print(t_detr)
print(Bmal_DetrData)

figure()
plot(t_detr, signal)

# DWT parameters
c_wavelet = "db20"      # wavelet
NumLevels = 7           # number of details considered

# Left/right margin handling

disregard_right = 12.
disregard_left  = disregard_right - 5.

## analyze first half of experiment
#disregard_left  = 12. - 5.
#disregard_right = t_detr[-1]-60.

## analyze second half of experiment
#disregard_left  = 60. - 5.
#disregard_right = t_detr[-1] - 120

# Apply MRA
LastOne = len(t_detr)
t_detr = array(t_detr).copy()[:LastOne]
level = pywt.swt_max_level(LastOne)

c_modwt = modwt(signal, c_wavelet, NumLevels)
c_wtmra = modwtmra(c_modwt, c_wavelet)

# Plot MRA results
PlotMRA_CoarseGrained_split(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, l_edge=disregard_left, r_edge=disregard_right)

print(Calculate_CoarseGrained_MRA_coefficients_split(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, l_edge=disregard_left, r_edge=disregard_right))

plt.savefig("./Results/Plots/Fig1C_GIMENExample.png")
plt.savefig("./Results/Plots/Fig1C_GIMENExample.svg")

### 
# Loop over all data
###

exp_index   = []
reporter    = []
c_cellline   = []

Noise       = []
Ultradian   = []
Circadian   = []
Infradian   = []

Noise_1st       = []
Ultradian_1st   = []
Circadian_1st   = []
Infradian_1st   = []

Noise_2nd       = []
Ultradian_2nd   = []
Circadian_2nd   = []
Infradian_2nd   = []

for ReporterDependentData, ReporterDependentCellines, c_rep in zip([Bmal_DetrData, Per2_DetrData], [Bmal_Celline, Per2_Celline], [BmalRep, Per2Rep]):
    t_detr = ReporterDependentData.index.to_numpy()
    #print(t_detr)
    for c_exp in ReporterDependentData.columns.to_list():
        exp_index.append(c_exp)
        str_split = c_exp.split("_")
        reporter.append(str_split[2])
        c_cellline.append(str_split[1])

        # downsample signal to obtain a nice circadian
        signal = ReporterDependentData[c_exp].to_numpy()[::3]

        # Apply MRA
        LastOne = len(t_detr)
        t_detr = array(t_detr).copy()[:LastOne]
        level = pywt.swt_max_level(LastOne)

        c_modwt = modwt(signal, c_wavelet, NumLevels)
        c_wtmra = modwtmra(c_modwt, c_wavelet)
        c_energ = Calculate_CoarseGrained_MRA_coefficients_split(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, l_edge=12.-5, r_edge=12.)
        c_energ_1st = Calculate_CoarseGrained_MRA_coefficients_split(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, l_edge=12.-5, r_edge=t_detr[-1]-60.)
        c_energ_2nd = Calculate_CoarseGrained_MRA_coefficients_split(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, l_edge=60. - 5., r_edge=t_detr[-1] - 120)
        #print(c_exp, c_exp.split("_"),  c_energ)
        Noise.append(c_energ[0])
        Ultradian.append(c_energ[1])
        Circadian.append(c_energ[2])
        Infradian.append(c_energ[3])

        Noise_1st.append(c_energ_1st[0])
        Ultradian_1st.append(c_energ_1st[1])
        Circadian_1st.append(c_energ_1st[2])
        Infradian_1st.append(c_energ_1st[3])

        Noise_2nd.append(c_energ_2nd[0])
        Ultradian_2nd.append(c_energ_2nd[1])
        Circadian_2nd.append(c_energ_2nd[2])
        Infradian_2nd.append(c_energ_2nd[3])

df = pd.DataFrame()
df["Cellline"]  = c_cellline
df["Reporter"]  = reporter
df["Noise"]     = Noise
df["Ultradian"] = Ultradian
df["Circadian"] = Circadian
df["Infradian"] = Infradian

df["Noise_1st"]     = Noise_1st
df["Ultradian_1st"] = Ultradian_1st
df["Circadian_1st"] = Circadian_1st
df["Infradian_1st"] = Infradian_1st

df["Noise_2nd"]     = Noise_2nd
df["Ultradian_2nd"] = Ultradian_2nd
df["Circadian_2nd"] = Circadian_2nd
df["Infradian_2nd"] = Infradian_2nd


df.index = exp_index
df.index.name='Experiment'

df.to_csv("./Results/Data/MRAAnalysis.csv", sep=",", index=True)

import seaborn as sns
g = sns.catplot(
    data=df, kind="bar",
    x="Cellline", y="Circadian", hue="Reporter",
    errorbar="sd", palette="dark", alpha=.6, height=6
)
#g.despine(left=True)
g.set_axis_labels("", "Circadianicity (%)")
g.legend.set_title("")
tight_layout()
#print(df)

savefig("./Results/Plots/MRA_BarPlotRhythmicity.svg", dpi=450)


figure(figsize=(6.4/1.5, 4.8/1.5))
# array with all cellines
tmp_celline = np.unique(df["Cellline"].values)
# dict that assigns color index to cellline
Color_Dict = {}
for i, j in enumerate(tmp_celline):
    Color_Dict[j] = i
print(Color_Dict)

# genrate data subsets based on reporter
dfB = df.loc[df["Reporter"]=="BMAL1"]
dfP = df.loc[df["Reporter"]=="PER2"]

cmap = matplotlib.cm.get_cmap('tab20')
print(cmap)

for n, i in enumerate(tmp_celline):
    if i in dfB["Cellline"].values and i in dfP["Cellline"].values:
        #print( i )
        tmp_B = dfB.loc[dfB["Cellline"] == i]
        tmp_P = dfP.loc[dfP["Cellline"] == i]
        errorbar(mean(tmp_B["Circadian"].values), mean(tmp_P["Circadian"].values), xerr=std(tmp_B["Circadian"].values), yerr=std(tmp_P["Circadian"].values), color=cmap(n), label=i, marker="o")
legend(loc=0, ncol=2)
xlabel("Bmal1 circadianicity (%)")
ylabel("Per2 circadianicity (%)")
plot([0, 100], [0, 100], "k--")
xlim(0, 100)
ylim(0, 100)

tight_layout()


savefig("./Results/Plots/MRA_Bmal1VSPer2.svg", dpi=450)



figure(figsize=(6.4/1.5, 4.8/1.5))
for n, i in enumerate(tmp_celline):
    tmp_B = dfB.loc[dfB["Cellline"] == i]
    tmp_P = dfP.loc[dfP["Cellline"] == i]
    errorbar(mean(tmp_B["Circadian_1st"].values), mean(tmp_B["Circadian_2nd"].values), xerr=std(tmp_B["Circadian_1st"].values), yerr=std(tmp_B["Circadian_2nd"].values), color=cmap(n), label=i, marker="o")
plot([0, 100], [0, 100], "k--")
xlabel("Bmal1 circadianicity 1st half(%)")
ylabel("Bmal1 circadianicity 2nd half(%)")
xlim(0, 100)
ylim(0, 100)
legend(loc=0, ncol=2)
tight_layout()

savefig("./Results/Plots/MRA_Bmal1_1stVS2nd.svg", dpi=450)


figure(figsize=(6.4/1.5, 4.8/1.5))
for n, i in enumerate(tmp_celline):
    tmp_B = dfB.loc[dfB["Cellline"] == i]
    tmp_P = dfP.loc[dfP["Cellline"] == i]
    errorbar(mean(tmp_P["Circadian_1st"].values), mean(tmp_P["Circadian_2nd"].values), xerr=std(tmp_P["Circadian_1st"].values), yerr=std(tmp_P["Circadian_2nd"].values), color=cmap(n), label=i, marker="o")
plot([0, 100], [0, 100], "k--")
xlabel("Per2 circadianicity 1st half(%)")
ylabel("Per2 circadianicity 2nd half(%)")
xlim(0, 100)
ylim(0, 100)
legend(loc=0, ncol=2)
tight_layout()
savefig("./Results/Plots/MRA_Per2_1stVS2nd.svg", dpi=450)

show()
