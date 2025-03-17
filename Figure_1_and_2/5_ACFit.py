
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

# import acf functions
import statsmodels.tsa.stattools as st


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


# Example analysis
c_data = Bmal_DetrData[Bmal_DetrData.columns[10]]

plot(t_detr, c_data)

figure()
acovf_data = st.acovf(c_data)
acf_data = st.acf(c_data, nlags=len(t_detr))

#plot(t_detr, acovf_data)
plot(t_detr-t_detr[0], acf_data)
#plot(acf_data)

# FIT DATA

#from lmfit import minimize, Parameters, fit_report

#def residual_damped(params, x, data=None):
    #D = params['D'].value
    #L = params['L'].value
    ##freq = params['freq'].value
    #tau = params['tau'].value
    #acovf_model = D/L*np.exp(-L*x)*np.cos(2.*np.pi/tau*x)
    #return (acovf_model - data)

#params = Parameters()
#params.add('D', value=0.02)
#params.add('L', value=0.05)
#params.add('tau', value=24.)

#print(sum(np.isnan(t_detr)), sum(np.isnan(acovf_data)),  sum(np.isnan(residual_damped(params, t_detr, acovf_data))))

##out = minimize(residual_damped, params, args=(t_detr, acovf_data), nan_policy='omit')
#out = minimize(residual_damped, params, args=(t_detr,), kws={'data': acovf_data})


###
# Fit using curve_fit
###

from scipy.optimize import curve_fit

def damped_func(x, A, L, tau):
    acovf_model = A*np.exp(-L*x)*np.cos(2.*np.pi/tau*x)
    return acovf_model

def damped_acovfunc(x, D, L, tau):
    acovf_model = D/L*np.exp(-L*x)*np.cos(2.*np.pi/tau*x)
    return acovf_model

print(damped_func(t_detr, 0.02, 0.05, 24.))
#plot(t_detr, damped_func(t_detr, 500., 0.01, 24.))

popt, pcov = curve_fit(damped_func, t_detr-t_detr[0], acf_data, p0=[500., 0.0005, 24.])
plot(t_detr-t_detr[0], damped_func(t_detr-t_detr[0], popt[0], popt[1], popt[2]), linestyle="--")

#params = Parameters()
#params.add('D', value=0.02)
#params.add('L', value=0.05)
#params.add('tau', value=24.)

#print(sum(np.isnan(t_detr)), sum(np.isnan(acovf_data)),  sum(np.isnan(residual_damped(params, t_detr, acovf_data))))

#out = minimize(residual_damped, params, args=(t_detr, acovf_data), nan_policy='omit')
#out = minimize(residual_damped, params, args=(t_detr,), kws={'data': acovf_data})
show()    

###
# ACF analysis
###
exp_index   = []
reporter    = []
c_cellline   = []

NoiseStrength   = []
DampingCoeff    = []
Period          = []
Amplitude       = []

for ReporterDependentData, ReporterDependentCellines, c_rep in zip([Bmal_DetrData, Per2_DetrData], [Bmal_Celline, Per2_Celline], [BmalRep, Per2Rep]):
    #figure(figsize=(6.4*2, 3.5*2))
    for i, cellline in enumerate(Celllines):
        #subplot(3,4,i+1)
        print(i, cellline)
        
        #if i in [0, 4, 8]:
            #ylabel(c_rep + " intensity")
        #if i in [8, 9, 10, 11]:
            #xlabel("time (h)")
        index_of_cellline = np.argwhere(ReporterDependentCellines == cellline).flatten()
        if len(index_of_cellline) != 0:
            Subset = ReporterDependentData.iloc[:, index_of_cellline]
            #max_power = 0
            for j, m in enumerate(Subset.columns.to_list()):
                print(m, m.split("_")[2], cellline)
                exp_index.append(m)
                reporter.append(m.split("_")[2])
                c_cellline.append(cellline)
                c_data = Subset[m]
                #acf_data = st.acf(c_data, nlags=len(t_detr))
                #popt, pcov = curve_fit(damped_func, t_detr-t_detr[0], acf_data, p0=[1., 0.005, 24.])
                
                acovf_data = st.acovf(c_data)
                L_ = 0.005
                #A_list.append( sqrt(out.params['D'].value/out.params['L'].value) )
                D_ = max(acovf_data)**2 * L_
                #print("Initial D: ", D_, max(acovf_data))
                popt, pcov = curve_fit(damped_acovfunc, t_detr-t_detr[0], acovf_data, p0=[D_, L_, 24.])
                NoiseStrength.append(popt[0])
                DampingCoeff.append(popt[1])
                Period.append(popt[2])
                Amplitude.append(np.sqrt(popt[0]/popt[1]))
                #wAn = pb.WAnalyzer(periods, dt, time_unit_label='hours')
                #modulus, transform = wAn.compute_spectrum(c_data, do_plot=False)
                #modulus, transform = wAn.compute_spectrum(c_data, do_plot=True)
                #ridge_results = wAn.get_maxRidge()
                #wAn.draw_Ridge()
                figure()
                #plot(t_detr-t_detr[0], acf_data, label="Data")
                #plot(t_detr-t_detr[0], damped_func(t_detr-t_detr[0], popt[0], popt[1], popt[2]), linestyle="--", label="Fit")
                plot(t_detr-t_detr[0], acovf_data, label="Data")
                plot(t_detr-t_detr[0], damped_acovfunc(t_detr-t_detr[0], popt[0], popt[1], popt[2]), linestyle="--", label="Fit")
                xlabel("Time lag (h)")
                #ylabel("Autocorrelation")
                ylabel("Autocovarianz")
                tight_layout()
                #savefig(f"./Results/Plots/ACF/Fits/{c_rep}_{m}.png")
                savefig(f"./Results/Plots/{c_rep}_{m}.png")
                
                #av_spec = wAn.get_averaged_spectrum()
                #print(ridge_results)
                #if c_rep == BmalRep:
                    #Bmal_Period[m]      = ridge_results["periods"].to_numpy()
                    #Bmal_Phase[m]       = ridge_results["phase"].to_numpy()
                    #Bmal_Amplitude[m]   = ridge_results["amplitude"].to_numpy()
                    #Bmal_Power[m]       = ridge_results["power"].to_numpy()
                    #Bmal_AverageSpec[m] = av_spec
                #elif c_rep == Per2Rep:
                    #Per2_Period[m]      = ridge_results["periods"].to_numpy()
                    #Per2_Phase[m]       = ridge_results["phase"].to_numpy()
                    #Per2_Amplitude[m]   = ridge_results["amplitude"].to_numpy()
                    #Per2_Power[m]       = ridge_results["power"].to_numpy()
                    #Per2_AverageSpec[m] = av_spec
                #print(Bmal_Period)
                #print(ridge_results)
                #if max(power_givenTLS)>max_power:
                    #max_power = max(power_givenTLS)
                #plot(t_detr, ridge_results["periods"], label=m)
                #plot(t, Subset[m], label=m)
            #legend(loc=0, prop={"size":5})
            #title(cellline)
            #xlim(0, t_detr[:-1])
            #ylim(0, max_power + 0.1 * max_power)
        else:
            #axis("off")
            pass


df = pd.DataFrame()
df["Cellline"]  = c_cellline
df["Reporter"]  = reporter
df.index = exp_index
df.index.name='Experiment'

df["NoiseStrength"] = NoiseStrength
df["DampingCoeff"] = DampingCoeff
df["Period"] = Period
df["Amplitude"] = Amplitude

df.to_csv("./Results/Data/ACovFAnalysis.csv", sep=",", index=True)



#import seaborn as sns

#g = sns.catplot(
    #data=df, kind="bar",
    #x="Cellline", y="Period", hue="Reporter",
    #errorbar="sd", palette="dark", alpha=.6, height=6
#)
##g.despine(left=True)
#g.set_axis_labels("", "Period (h)")
#g.legend.set_title("")
#tight_layout()
##print(df)

#savefig("./Results/Plots/ACovF/Periods.png", dpi=450)


#g = sns.catplot(
    #data=df, kind="bar",
    #x="Cellline", y="DampingCoeff", hue="Reporter",
    #errorbar="sd", palette="dark", alpha=.6, height=6
#)
##g.despine(left=True)
#g.set_axis_labels("", "Period (h)")
#g.legend.set_title("")
#tight_layout()
##print(df)

#savefig("./Results/Plots/ACovF/DampingCoefficient.png", dpi=450)
