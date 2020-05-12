#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:36:28 2020

@author: s.bykov
"""

import os
import seaborn as sns
from PipelineNuSTAR.core import *


#%% simulate lc #https://stingray.readthedocs.io/en/latest/notebooks/Simulator/Simulator%20Tutorial.html
import stingray
from stingray import Lightcurve, Crossspectrum, Powerspectrum
from stingray.simulator import simulator
from Misc.TimeSeries import cross_correlation

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_0.1')



def simulate_lc712_from_powerspectrum(N=10000,mean=258.2,rms=0.06,dt=0.1,
                                   plot_results=0):

    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    #from xspec
    xspec_pars=[0.743175,
 0.0937313,
 0.0303734,
  0.182193,
 0.0725758,
  0.323264,
  0.228514,
     0.001,
 0.0236255,
 0.0559709,
 0.0245674,
  0.203363,
         0,
  0.716638,
    1.1058,
 0.0426574]

    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node191.html
    def lorentzian( x, x0, gam, a ):
        #return a * gam**2 / ( gam**2 + ( x - x0 )**2)
        return a*(gam/(2*np.pi))/( (x-x0)**2 + (gam/2)**2  )
    def gaussian(x,x0,sigma,N):
        return N*1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-x0)**2/(2*sigma**2))
    def po(x,gamma,N):
        return N*x**(-gamma)
#model  gaussian + gaussian + lorentz + gaussian + powerlaw + powerlaw

    spectrum = gaussian(w, *xspec_pars[0:3])+gaussian(w, *xspec_pars[3:6])+ lorentzian(w, *xspec_pars[6:9])+gaussian(w,*xspec_pars[9:12])+po(w,*xspec_pars[12:14])+po(w,*xspec_pars[14:16])


    lc = sim.simulate(spectrum)
    #tmp=np.random.normal(lc.counts-lc.counts,270)
    #lc.counts=lc.counts+np.random.normal(lc.counts-lc.counts,270)
    if plot_results:
        plt.figure()
        plt.plot(lc.time,lc.counts)
        plt.show()

        plt.figure()
        ps=Powerspectrum(lc)
        ps=ps.rebin_log(0.01)
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

ccfs=[]
for i in range(500):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    lc712=lc712.rebin(10)
    A=0.8
    deltaT=10
    lc67=iron_band_model_lc(lc712, A, deltaT)
    ccf=find_ccf(lc712, lc67,deltaT=deltaT,A=A,plot=0)
    ccfs.append(ccf.ccf)
ccfs=np.asarray(ccfs)


#%%plot simulations
fig,ax_ccf=plt.subplots(figsize=(16,6))

ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0),ccfs.std(axis=0),label=f'simulations dT={deltaT}s; A={A}')
ax_ccf.set_xlim(-75,75)
ax_ccf.set_ylim(-0.15,1)


def plot_ccf(filepath,ax):
    ccf=np.genfromtxt(filepath,skip_header=3)
    N=int(ccf[:,0].shape[0]/2)
    ax.errorbar(ccf[:,0],ccf[:,2],ccf[:,3],drawstyle='steps-mid',label='data')

    #ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')
    ax.errorbar(ccf[N:,0],ccf[0:N+1:,2][::-1],ccf[0:N+1:,3][::-1],alpha=0.5,color='m',drawstyle='steps-mid',label='data (neg delay)')
    ax.set_xlabel('Delay, s')

    ax.set_ylabel('CCF')
    fig.tight_layout()
    sns.despine(fig,top=1,right=0)
    plt.show()

plot_ccf('../lc712/ccf_712vs67.qdp',ax_ccf)

ax_ccf.set_xlabel('iron delay, s')
ax_ccf.legend()
plt.show()

plt.savefig(f'simulations/ccf_dt{deltaT}s_A{A}.png')
