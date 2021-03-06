#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 16:02:07 2020

@author: s.bykov
"""

#%% imports
from PipelineNuSTAR.core import *
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

path_to_lc='/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products'

class TimeSeries():

    def __init__(self,lcname):
        self.fits=fits.open(path_to_lc+f'/{lcname}/{lcname}AB_sr.lc_bary')
        self.time=self.fits[1].data['time']
        self.rate=self.fits[1].data['rate']
        self.error=self.fits[1].data['error']
        self.binsize=np.median(np.diff(self.time))
        self.fits.close()
    def divide(self,val):
        self.rate=self.rate/val
        self.error=self.error/val

stop
#%% plot arf
arf=fits.open('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/spe_and_lc/spe_and_lcA_sr.arf')
spe=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc46/spe_data.qdp',skip_header=3)

enaxis=arf[1].data['ENERG_LO']
rsp=arf[1].data['SPECRESP']
plt.plot(enaxis,rsp,label=' arf')
plt.xlim(4,9)

en=np.array([5,6.5,8])
enerr=np.array([1,0.5,1])
dE=np.array([2,1,2])


#for mean_e in en:
#    plt.axvline(mean_e,color='k')

#plt.show()

#area=np.array([230,310,347])

#plot(enaxis,rsp*cutoffpl(enaxis,norm=0.170091*(51/36)),label='cutoff pl orig model from xspec * arf, scaled')

#%%read lcurves
lc46=TimeSeries('lc46')
lc67=TimeSeries('lc67')
lc79=TimeSeries('lc79')
lc712=TimeSeries('lc712')


for k,lc in enumerate([lc46,lc67,lc79]):
    lc=lc.divide(dE[k])

rate=np.array([lc.rate.mean() for lc in [lc46,lc67,lc79]])
rate_error=np.array([lc.rate.std() for lc in [lc46,lc67,lc79]])

#%% plot spectra etc
plt.figure()
# def cutoffpl(x,gamma=-0.337940,efold=7.91525,norm=0.170091):
#     return norm*x**(-gamma)*np.exp(-x/efold)
# plot(enaxis,rsp*cutoffpl(enaxis,norm=0.170091*(51/36)),label='cutoff pl orig model from xspec * arf, scaled')

a,b=np.polyfit(en[[0,2]],rate[[0,2]],1)
#enaxis=np.linspace(4,9,100)
plt.plot(enaxis,a*enaxis+b,label='linear best model ignoring 6-7 kev')

plt.errorbar(en,rate,rate_error,dE/2,label='spectra_From_lc')
plt.plot(spe[:,0],spe[:,2]*63/46,'k.',label='data from xspec + 36%')
#plt.xscale('log')
plt.show()
plt.ylabel(' counts/s / keV')
plt.xlabel('energy')
plt.legend()
plt.xlim(4,9)
plt.ylim(30,80)

def best_line(x,N):
    return x*a+N


arf.close()



#%% plot a few points
N=len(lc46.time)


fig,ax=plt.subplots(3)
from random import randint
randi=[randint(0, N) for p in range(0, 3)]
for k,i in enumerate(randi):
    rate=np.array([lc46.rate[i],lc67.rate[i],lc79.rate[i]])
    rate_error=np.array([lc46.error[i],lc67.error[i],lc79.error[i]])

    N_opt,N_opt_err=curve_fit(best_line,en[[0,2]],rate[[0,2]],
                              p0=b,sigma=rate_error[[0,2]],absolute_sigma=True)
    print(N_opt_err)
    N_opt_err=np.sqrt(N_opt_err[0])


    myline=lambda x: a*x+N_opt
    ax[k].plot(enaxis,myline(enaxis),'k:',alpha=0.3)

    ax[k].errorbar(en,rate,rate_error,dE/2,label='spectra_From_lc')


    diff=rate[1]-myline(en[1])
    diff_err=np.sqrt(rate_error[1]**2+N_opt_err**2)
    #diff_err=np.sqrt(rate_error[1]**2)

    ax[k].vlines(en[1],rate[1],rate[1]-diff,color='r',alpha=0.7)

    print(i,diff,diff_err)



plt.show()
for axx in ax:
    axx.set_ylabel(' counts/s / keV')
    axx.set_xlabel('energy')
    #plt.legend()
    axx.set_xlim(4,9)
    axx.set_ylim(30,80)



#%% find fe flux in time

def find_fe_flux(lc46,lc67,lc79):
    N=len(lc46.rate)

    diff=[]
    diff_err=[]
    for i in range(N):
        f=np.array([lc46.rate[i],lc67.rate[i],lc79.rate[i]])
        ferr=np.array([lc46.error[i],lc67.error[i],lc79.error[i]])


        N_opt,N_opt_err=curve_fit(best_line,en[[0,2]],f[[0,2]],
                              p0=b,sigma=ferr[[0,2]],absolute_sigma=True)
        N_opt_err=np.sqrt(N_opt_err[0])

        myline=lambda x: a*x+N_opt

        fe_flux=f[1]-myline(en[1])
        fe_flux_err=np.sqrt(ferr[1]**2+N_opt_err**2)
        #fe_flux_err=np.sqrt(ferr[1]**2)

        diff.append(fe_flux)
        diff_err.append(fe_flux_err)

    diff=np.asarray(diff)
    diff=np.reshape(diff,N)
    diff_err=np.asarray(diff_err)
    diff_err=np.reshape(diff_err,N)
    return diff,diff_err

fe_flux,fe_flux_err=find_fe_flux(lc46, lc67, lc79)

Flux_tot=lc67.rate
Flux_fe=fe_flux
Flux_cont=Flux_tot-Flux_fe

eqw=(Flux_tot-Flux_cont)*dE[1]/Flux_cont*1000


#%% plot iron intensity
figure()
plt.errorbar(lc46.time,fe_flux,fe_flux_err)
plt.xlabel('Time, s')
plt.ylabel('Iron line flux')
plt.show()

figure()
plt.hist(fe_flux,bins=50,label='Fe_flux')
plt.xlabel('Iron line flux ')
plt.show()

print(f'Mean significance: {np.mean(fe_flux/fe_flux_err)}')
print(f'Mean Flux: {np.mean(fe_flux)}')
print(f'Std Flux: {np.std(fe_flux)}')
print(f'Mean Error: {np.mean(fe_flux_err)}')

'''
with additional error
Mean significance: 1.8413464247438123
Mean Flux: 9.163717277556117
Std Flux: 3.504459336753297
Mean Error: 4.9813509551603365

without additional error
Mean significance: 2.0250680158870975
Mean Flux: 9.163717277556117
Std Flux: 3.504459336753297
Mean Error: 4.5219418245787715

'''


#%% iron line flux save in fits

STOP
with fits.open('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/fe_line_sr.lc_bary', mode='update') as hdul:
    hdul[1].data['rate']=fe_flux
    hdul[1].data['error']=fe_flux_err
    hdul.flush()  # changes are written back to original.fits

with fits.open('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/fe_line_eqw.lc_bary', mode='update') as hdul:
    hdul[1].data['rate']=eqw
    hdul[1].data['error']=eqw/2
    hdul.flush()  # changes are written back to original.fits


#%% read cross corr output
matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.9
plt.subplots_adjust(wspace=2)
plt.subplots_adjust(hspace=1)


fig,ax=plt.subplots()
ccf=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/ccf.qdp',skip_header=3)
N=int(ccf[:,0].shape[0]/2)
ax.errorbar(ccf[:,0],ccf[:,2],ccf[:,3],drawstyle='steps-mid')

ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')

#ax.set_xlabel('Iron line leads <---Delay, s--->  7-12 keV Flux leads')
ax.set_xlabel('Delay, s')

ax.set_ylabel('CCF \n Iron line flux vs FLux 7-12 keV')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()



#%% xspec fluxes for mean observation
'''
XSPEC12>flux 6 7
 Model Flux   0.14894 photons (1.5498e-09 ergs/cm^2/s) range (6.0000 - 7.0000 keV)

Model cflux<1>*gaussian<2> +  cutoffpl<3> Source No.: 1   Active/On
Model Model Component  Parameter  Unit     Value
 par  comp
   1    1   cflux      Emin       keV      3.00000      frozen
   2    1   cflux      Emax       keV      12.0000      frozen
   3    1   cflux      lg10Flux   cgs      -10.0240     +/-  1.44138E-02
   4    2   gaussian   LineE      keV      6.51042      +/-  1.44220E-02
   5    2   gaussian   Sigma      keV      0.300000     frozen
   6    2   gaussian   norm                9.07066E-03  frozen
   7    3   cutoffpl   PhoIndex            -0.351107    +/-  2.58053E-03
   8    3   cutoffpl   HighECut   keV      7.84568      +/-  3.85120E-02
   9    3   cutoffpl   norm                0.167056     frozen
____________________________________________________________________

10**-10.0240
Out[11]: 9.462371613657949e-11

10**-10.0240/1.5498e-09
Out[12]: 0.06105543691868596

'''


#%% simulation
STOP
#%% simulation of sine waves
T=4000
binsize=10
N=int(T/binsize)
time=np.linspace(0,T,N) #typical nustar window duration of 4000 sec with 10 sec binsize

###y1 - 7-12 keV
P=150
y1_mean=250
y1_error=20
y1_ampt=50

y1=y1_mean+y1_ampt*np.sin(time*2*np.pi/P)
y1+=np.random.normal(time-time,y1_error)
y1_err=y1_error


###y2 - iron line
deltaT=10
y2_mean=10
y2_error=5
y2_ampt=3

y2=y2_mean+y2_ampt*np.sin((time-deltaT)*2*np.pi/P)
y2+=np.random.normal(time-time,y2_error)
y2_err=y2_error

#%% plots
fig,axs=plt.subplots(3,)#sharex='all')
axs[0].errorbar(time,y1,y1_err)
axs[1].errorbar(time,y2,y2_err)
ax=axs[2]
#plt.show()



#%%crosscorr
from Misc.TimeSeries.cross_correlation import cross_correlation


lag,corr=cross_correlation(time, y1, y2,circular=0)
N=int(lag.shape[0]/2)


#fig,ax=plt.subplots()

ax.plot(lag,corr,label='CCF \n y1=7-12 keV \n y2=Iron flux')

ax.plot(-lag[N:],corr[N:],alpha=0.5,color='r')



ax.set_xlim(-150,150)
ax.set_xlabel('y2 lags <--0--> y1 lags')
ax.set_ylabel('CCF')
ax.legend()

axs[0].set_title(f'''
             y2 deltaT={deltaT}
             y2_mean={y2_mean}; y2_error={y2_error};y2_ampitude={y2_ampt}
             ''',fontsize=12)
plt.show()

savepath='/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/simulations/'
plt.savefig(savepath+f'ccf_{deltaT}.png')





#%% simulation interpolate


lcdata=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/data_lc712_segm4.qdp',skip_header=0)
time=lcdata[:,0]
time-=time[0]
rate=lcdata[:,2]
error=lcdata[:,3]
from scipy import interpolate
from scipy.signal import savgol_filter


itp = interpolate.interp1d(time,rate, kind='linear')
window_size, poly_order =21, 3
yy_sg = savgol_filter(itp(time), window_size, poly_order)


#f = interpolate.interp1d(time, rate, kind='cubic')
#f=interpolate.CubicSpline(time,rate)
plt.errorbar(time,rate,error,alpha=0.6)
#timeaxis=np.linspace(time[0],time[-1],2000)
plt.plot(time, yy_sg, 'r-', label= "Smoothed curve")
#plt.plot(timeaxis,f(timeaxis),'r-',lw=3)

f=interpolate.interp1d(time, yy_sg,fill_value='extrapolate')

plt.show()

#%% create data

y1=rate
y1_err=error


deltaT=20
y2_mean=10
y2_error=5
y2_ampl=10

y2=f((time-deltaT))/mean(y1)*y2_mean*y2_ampl
y2+=np.random.normal(time-time,y2_error)
y2_err=y2_error


#%% plots
fig,axs=plt.subplots(3,)#sharex='all')
axs[0].errorbar(time,y1,y1_err)
axs[0].plot(time,f((time)),'r:',zorder=10)

axs[1].errorbar(time,y2,y2_err)
axs[1].plot(time,f((time-deltaT))/mean(y1)*y2_mean*y2_ampl,'r:',zorder=10)
ax=axs[2]
#plt.show()



#%%crosscorr
from Misc.TimeSeries.cross_correlation import cross_correlation


lag,corr=cross_correlation(time, y1, y2,circular=0)
N=int(lag.shape[0]/2)


#fig,ax=plt.subplots()

ax.plot(lag,corr,label='CCF \n y1=7-12 keV \n y2=Iron flux')

ax.plot(-lag[N:],corr[N:],alpha=0.5,color='r')



ax.set_xlim(-150,150)
ax.set_xlabel('y2 lags <--0--> y1 lags')
ax.set_ylabel('CCF')
ax.legend()

axs[0].set_title(f'''
             y2 deltaT={deltaT}
             y2_mean={y2_mean}; y2_error={y2_error};y2_ampitude={y2_ampl}
             ''',fontsize=12)
plt.show()

savepath='/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712/simulations/'
plt.savefig(savepath+f'ccf_{deltaT}s_{y2_ampl}ampl_real_lc.png')
