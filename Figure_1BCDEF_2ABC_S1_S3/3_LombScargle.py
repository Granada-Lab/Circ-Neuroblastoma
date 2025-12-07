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


# define reporter names
BmalRep = "Bmal1-Luc"
Per2Rep = "Per2::LUC"


# set up detrending parameters
dt          = 10./60     # in [h]
T_cut_off   = 48  # cut-off period for sinc-filter
skip        = int(5 // dt)



#Bmal_DetrData.to_csv("./Results/Detrended_Bmal.csv", sep=";")
#Per2_DetrData.to_csv("./Results/Detrended_Per2.csv", sep=";")

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
# Lomb-Scargle analysis
###

num_rows    = 4
num_columns = 4


TLS = arange(10., 50, 0.1)
for ReporterDependentData, ReporterDependentCellines, c_rep in zip([Bmal_DetrData, Per2_DetrData], [Bmal_Celline, Per2_Celline], [BmalRep, Per2Rep]):
    figure(figsize=(6.4*2, 3.5*2))
    for i, cellline in enumerate(Celllines):
        subplot(num_rows, num_columns,i+1)
        print(i, cellline)
        #if i in [0, 4, 8]:
            #ylabel(c_rep + " intensity")
        #if i in [8, 9, 10, 11]:
            #xlabel("time (h)")
        index_of_cellline = np.argwhere(ReporterDependentCellines == cellline).flatten()
        if len(index_of_cellline) != 0:
            Subset = ReporterDependentData.iloc[:, index_of_cellline]
            max_power = 0
            for j, m in enumerate(Subset.columns.to_list()):
                c_data = Subset[m]
                ls = LombScargle(t_detr, c_data)
                power_givenTLS    = ls.power(1./TLS)
                if max(power_givenTLS)>max_power:
                    max_power = max(power_givenTLS)
                plot(TLS, power_givenTLS, label=m)
                #plot(t, Subset[m], label=m)
            legend(loc=0, prop={"size":5})
            title(cellline)
            xlim(0, TLS[-1])
            ylim(0, max_power + 0.1 * max_power)
        else:
            axis("off")
    tight_layout()
    #savefig("Raw"+c_rep+".png", dpi=450)



####
## Multiresolution analysis (MRA)
####

#def PlotMRA_CoarseGrained(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=12.):
    #f, axarr = plt.subplots(5, sharex=True, figsize=(6/1.7, 12/2.5))
    #c_max = ceil(max(signal.flatten()))

    #i_edge = round(t_edge/dt)

    #OverallEnergy = sum([sum(n[i_edge:-i_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-i_edge])
    #D1Energy = sum(MRA[0][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D2Energy = sum(MRA[1][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D3Energy = sum(MRA[2][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D4Energy = sum(MRA[3][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D5Energy = sum(MRA[4][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D6Energy = sum(MRA[5][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D7Energy = sum(MRA[6][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    #D8Energy = sum(MRA[7][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy


    #D12  = MRA[0] + MRA[1]
    #D34  = MRA[2] + MRA[3]
    #D5   = MRA[4]
    #D678 = MRA[5] + MRA[6] + MRA[7]

    #c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    #c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy+D8Energy]
    #for i, j in enumerate([D12, D34, D5, D678]):
        #axarr[i].plot(t, j, color="gray", linestyle="--")

        #c_maxtmp = max(abs(j[i_edge:-i_edge]))
        #axarr[i].plot(t[i_edge:-i_edge], j[i_edge:-i_edge], label = str(round(c_energy[i]*100,2)))
        #axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        #axarr[i].set_xlim(0, t[-1])
        #axarr[i].set_ylim(-c_maxtmp, c_maxtmp)
        #axarr[i].set_ylabel(c_label[i])
        #leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        #leg2.get_frame().set_alpha(0.5)


    #axarr[4].plot(t, signal, "k-", linewidth=2, label="Signal")
    #axarr[4].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    #axarr[4].set_ylim(-c_max, c_max)
    #axarr[4].set_yticks([-c_max, 0, c_max])
    #leg1 = axarr[4].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    #leg1.get_frame().set_alpha(0.5)

    #axarr[4].set_xlabel("Time in h")
    #f.tight_layout()

## Downsample data to obtain a proper circadian frequency band
#t_detr = t_detr[::3]
#signal = Bmal_DetrData["20210309_SKNAS_BMAL1_3"].to_numpy()[::3]
#print(t_detr)
#print(Bmal_DetrData)

#figure()
#plot(t_detr, signal)

## DWT parameters
#c_wavelet = "db20"      # wavelet
#NumLevels = 7           # number of details considered

## Apply MRA
#LastOne = len(t)
#t_detr = array(t_detr).copy()[:LastOne]
#level = pywt.swt_max_level(LastOne)

#c_modwt = modwt(signal, c_wavelet, NumLevels)
#c_wtmra = modwtmra(c_modwt, c_wavelet)

## Plot MRA results
#PlotMRA_CoarseGrained(t_detr, signal, c_wtmra, dt=30./60, NumLevels=7, t_edge=12.)

show()    
