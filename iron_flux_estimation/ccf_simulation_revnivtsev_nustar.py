#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:36:28 2020

@author: s.bykov
"""


#%% simulate lc #https://stingray.readthedocs.io/en/latest/notebooks/Simulator/Simulator%20Tutorial.html
import os
import numpy as np
import seaborn as sns
from PipelineNuSTAR.core import *


import stingray
from stingray import Lightcurve, Crossspectrum, Powerspectrum
from stingray.simulator import simulator
from Misc.TimeSeries import cross_correlation

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712')



def simulate_lc712_from_powerspectrum(mean=253.4,rms=np.sqrt(402.5)/253,dt=10,
                                   plot_results=0):
    N=int(1000/dt)
    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    #from xspec
    xspec_pars=[
 0.00912596,
 0.00300072,
  0.0785102,
          0,
 0.00215227,
   0.240521,
          0,
    5.85424,]

    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node191.html
    def lorentzian( x, x0, gam, a ):
        #return a * gam**2 / ( gam**2 + ( x - x0 )**2)
        return a*(gam/(2*np.pi))/( (x-x0)**2 + (gam/2)**2  )
    def gaussian(x,x0,sigma,N):
        return N*1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-x0)**2/(2*sigma**2))
    def po(x,gamma,N):
        return N*x**(-gamma)
#model  const+lore+gauss

    spectrum = lorentzian(w, *xspec_pars[0:3])+lorentzian(w, *xspec_pars[3:6])+ po(w,*xspec_pars[6:8])


    lc = sim.simulate(spectrum)

    #tmp=np.random.normal(lc.counts-lc.counts,10)
    #lc.counts=lc.counts+tmp


    if plot_results:
        plt.figure()
        plt.plot(lc.time,lc.counts)
        plt.show()

        plt.figure()
        ps=Powerspectrum(lc)
        ps=ps.rebin_log(0.05)
        plt.loglog(ps.freq,ps.power)
        plt.show()
    return lc


def iron_band_model_lc(lc712,A,deltaT):
    #lc67(t)= A*lc712(t)+(1-A)*lc712(t+dt)
    #dt>0 : delay with respect to lc712
    lc712_cts=lc712.counts

    tmpcounts=np.roll(lc712_cts,int(deltaT/lc712.dt))
    lc67_cts=A*lc712_cts+(1-A)*tmpcounts
    lc67=Lightcurve(lc712.time, lc67_cts)
    return lc67

def find_ccf(lc712,lc67,plot=0,deltaT=0,A=0):
    CCF_obj_crosscorr=cross_correlation.CrossCorrelation(lc712.time, lc67.counts,lc712.counts,circular=0)
    CCF_obj_crosscorr.calc_ccf()
    #plt.close()

    if plot:
        fig,ax=plt.subplots(1,sharex='all')
        ax.plot(CCF_obj_crosscorr.lag,CCF_obj_crosscorr.ccf,'b:.',label=f'simulations dT={deltaT}s; A={A}')
        ax.grid()
        N=int(CCF_obj_crosscorr.lag.shape[0]/2)
        ax.set_xlabel('Delay, s')
        ax.set_ylabel('CCF')
        ax.legend()
        plt.xlim(-15,15)
        plt.ylim(-0.5,1)
        plt.show()
    return CCF_obj_crosscorr

#%% test simulator
simulate_lc712_from_powerspectrum(plot_results=1)


#%% simulate
ccfs=[]
for i in range(100):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    A=0.7
    deltaT=30
    lc67=iron_band_model_lc(lc712, A, deltaT)

    #lc712=lc712.rebin(1)
    #lc67=lc712.rebin(1)

    ccf=find_ccf(lc712, lc67,deltaT=deltaT,A=A,plot=0)
    ccfs.append(ccf.ccf)
ccfs=np.asarray(ccfs)


#%%plot simulations

matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.9
plt.subplots_adjust(wspace=2)
plt.subplots_adjust(hspace=1)

fig,ax_ccf=plt.subplots()

ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0)*0.25,ccfs.std(axis=0)*0.25,label=f'simulations dT={deltaT}s; A={A}',marker='.',color='k',alpha=0.7)
ax_ccf.set_xlim(-150,150)
#ax_ccf.set_ylim(-0.2,1.1)

def plot_data_ccf(name,ax_ccf):
    ccf_data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/{name}.qdp',skip_header=3)
    N=int(ccf_data[:,0].shape[0]/2)
    norm=1#np.max(ccf_data[:,2])
    ax_ccf.errorbar(ccf_data[:,0],ccf_data[:,2]/norm,ccf_data[:,3]/norm,drawstyle='steps-mid',alpha=0.8,label=name,color='b')
    #ax_ccf.errorbar(-ccf_data[:N+1,0],ccf_data[:N+1,2]/norm,ccf_data[:N+1,3]/norm,alpha=0.5,color='r',drawstyle='steps-mid')
    return None
#plot_data_ccf('ccf_0.8', ax_ccf)
plot_data_ccf('ccf_0.85', ax_ccf)
#plot_data_ccf('ccf_1', ax_ccf)

ax_ccf.set_xlabel('Delay, s')
ax_ccf.set_ylabel('CCF')
#ax_ccf.legend()
plt.show()

#plt.savefig(f'simulations/ccf_A{A}_dt{deltaT}s_powspec_simulations.pdf')
plt.savefig(f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/fe_line_intensity_0.85.pdf')