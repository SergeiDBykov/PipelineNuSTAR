#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:36:28 2020

@author: s.bykov
"""

stop
#%% read cross corr output
matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.9
plt.subplots_adjust(wspace=2)
plt.subplots_adjust(hspace=1)


fig,ax=plt.subplots()
ccf=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_0.1/ccf_orb_corr.qdp',skip_header=3)
N=int(ccf[:,0].shape[0]/2)
ax.errorbar(ccf[:,0],ccf[:,2],ccf[:,3],drawstyle='steps-mid')

ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')

#ax.set_xlabel('Iron line leads <---Delay, s--->  7-12 keV Flux leads')
ax.set_xlabel('Delay, s')

ax.set_ylabel('CCF \n Iron line flux vs FLux 7-12 keV')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()




#%% compare pdf
figure()
for filename in ['pds.qdp','../lc67_0.1/pds.qdp']:
    data=np.genfromtxt(filename,skip_header=3)
    plt.errorbar(data[:,0],data[:,2],data[:,3],data[:,1])

plt.xscale('log')
plt.yscale('log')
plt.show()

#%% stingray cross cpectra

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_0.1/')

lc1=fits.open('lc712AB_sr.lc_bary_orb_corr')
lc2=fits.open('../lc67_0.1/lc67AB_sr.lc_bary_orb_corr')
#lc2=fits.open('ch1819.lc')

from stingray import Lightcurve
lc1=Lightcurve(lc1[1].data['time'], lc1[1].data['rate'],
              err=lc1[1].data['error'],input_counts=0)
#lc1=lc1.rebin(0.15)

from stingray.events import EventList
#evt1=EventList.from_lc(lc1)
#evt1.simulate_energies([[5.5],[1]])

lc2=Lightcurve(lc2[1].data['time'], lc2[1].data['rate'],
              err=lc2[1].data['error'],input_counts=0)
#lc2=lc2.rebin(0.15)

#evt2=EventList.from_lc(lc2)
#evt2.simulate_energies([6.5,1])
#evt=evt1.join(evt2)

#%%pds
from stingray import Powerspectrum, AveragedPowerspectrum

avg_ps=AveragedPowerspectrum(lc1,100)


#%% cross
from stingray import Crossspectrum,AveragedCrossspectrum

avg_cs=AveragedCrossspectrum(lc1,lc2,20)

freq_lags, freq_lags_err = avg_cs.time_lag()


fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.hlines(0, avg_cs.freq[0], avg_cs.freq[-1], color='black', linestyle='dashed', lw=2)
ax.errorbar(avg_cs.freq, freq_lags, yerr=freq_lags_err,fmt="o", lw=1, color='blue')
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Time lag (s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(which='major', width=1.5, length=7)
ax.tick_params(which='minor', width=1.5, length=4)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)
plt.show()



#%% simulate lc #https://stingray.readthedocs.io/en/latest/notebooks/Simulator/Simulator%20Tutorial.html

from stingray import Lightcurve, Crossspectrum, sampledata
from stingray.simulator import simulator, models

#%% simulate 7-12 keV lightcurve with 4s binsize
sim = simulator.Simulator(N=4000, mean=200, dt=1)

figure()
w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]

def lorentzian( x, x0, gam, a ):
    return a * gam**2 / ( gam**2 + ( x - x0 )**2)


#spectrum = np.power((1/w),0.75)*0.26+1.27
spectrum=   np.power((1/w),0.98)*0.06+1.65 + lorentzian(w, 0.056710493, 0.027166799, 3.52411) +lorentzian(w, 0.009623506, 0.0029425912, 13.037277) + lorentzian(w,  0.21675608, 0.04158771 ,1.589633)
plt.loglog(w,spectrum)

lc = sim.simulate(spectrum)
time=lc.time

figure()
plt.plot(lc.counts)
plt.show()

#%% model parameters
# lc67(t)= (1-A)*lc712(t)+A*lc712(t+dt)
#dt>0 : delay with respect to lc712
lc712_cts=lc.counts


A=0.8
dt=2

tmpcounts=np.roll(lc712_cts,int(dt/lc712.dt))


lc67_cts=A*lc712_cts+(1-A)*tmpcounts


lc712=Lightcurve(time,lc712_cts)
lc67=Lightcurve(time, lc67_cts)

figure()
plt.plot(lc712.counts,'b.:')
plt.plot(lc67.counts,'rs-.')

plt.show()

#%% ccf stuff

from Misc.TimeSeries import cross_correlation


CCF_obj=cross_correlation.CrossCorrelation(time, lc712.counts, lc67.counts,circular=0)
CCF_obj.calc_ccf(subtract_mean=1)

plt.plot(CCF_obj.lag,CCF_obj.ccf,'g-',ms=1)
N=int(CCF_obj.lag.shape[0]/2)

plt.plot(-CCF_obj.lag[N:],CCF_obj.ccf[N:],alpha=0.5,color='r')


plt.xlim(-70,70)


