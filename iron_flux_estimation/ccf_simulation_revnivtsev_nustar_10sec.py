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

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712')



def simulate_lc712_from_powerspectrum(N=10000,mean=253.4,rms=0.08,dt=10,
                                   plot_results=0):

    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    #from xspec
    xspec_pars=[4.250518,
                0., +0.0030077815, 91.51196,
                0.009163746, 0.0028483302, 14.861738]

    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node191.html
    def lorentzian( x, x0, gam, a ):
        #return a * gam**2 / ( gam**2 + ( x - x0 )**2)
        return a*(gam/(2*np.pi))/( (x-x0)**2 + (gam/2)**2  )
    def gaussian(x,x0,sigma,N):
        return N*1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-x0)**2/(2*sigma**2))
    def po(x,gamma,N):
        return N*x**(-gamma)
#model  const+lore+gauss

    spectrum = xspec_pars[0]+lorentzian(w, *xspec_pars[1:4])+lorentzian(w, *xspec_pars[4:7])


    lc = sim.simulate(spectrum)
    tmp=np.random.normal(lc.counts-lc.counts,10)


    lc.counts=lc.counts+tmp
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
for i in range(50):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    A=0.8
    deltaT=10
    lc67=iron_band_model_lc(lc712, A, deltaT)
    ccf=find_ccf(lc712, lc67,deltaT=deltaT,A=A,plot=0)
    ccfs.append(ccf.ccf)
ccfs=np.asarray(ccfs)


#%%plot simulations
fig,ax_ccf=plt.subplots(figsize=(16,6))

ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0),ccfs.std(axis=0),label=f'simulations dT={deltaT}s; A={A}',marker='s')
ax_ccf.set_xlim(-75,75)
ax_ccf.set_ylim(-0.2,1.1)


ccf_data=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/ccf.qdp',skip_header=3)

ax_ccf.errorbar(ccf_data[:,0],ccf_data[:,2]/np.max(ccf_data[:,2]),ccf_data[:,3]/np.max(ccf_data[:,2]),drawstyle='steps-mid',label='data (normed)')


ax_ccf.set_xlabel('iron delay, s')
ax_ccf.legend()
plt.show()

#plt.savefig(f'simulations/ccf_dt{deltaT}s_A{A}_powspec_simulations.png')
