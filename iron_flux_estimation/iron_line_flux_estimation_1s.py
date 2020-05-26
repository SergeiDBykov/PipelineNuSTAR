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
        self.fits=fits.open(path_to_lc+f'/{lcname}_1/{lcname}AB_sr.lc_bary_orb_corr')
        self.time=self.fits[1].data['time']
        self.rate=self.fits[1].data['rate']
        self.error=self.fits[1].data['error']
        self.binsize=np.median(np.diff(self.time))
        self.fits.close()
    def divide(self,val):
        self.rate=self.rate/val
        self.error=self.error/val

stop

#%%read lcurves
en=np.array([5,6.5,8])
enerr=np.array([1,0.5,1])
dE=np.array([2,1,2])


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

spe=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc46/spe_data.qdp',skip_header=3)

a,b=np.polyfit(en[[0,2]],rate[[0,2]],1)
#enaxis=np.linspace(4,9,100)
enaxis=np.linspace(4,9,100)
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
    axx.set_ylim(30,120)



#%% find fe flux in time

def find_fe_flux(lc46,lc67,lc79,frac=1):
    N=len(lc46.rate)

    diff=[]
    diff_err=[]
    test_ratio=[]
    print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')
    for i in range(N):
        f=np.array([lc46.rate[i],lc67.rate[i],lc79.rate[i]])
        ferr=np.array([lc46.error[i],lc67.error[i],lc79.error[i]])


        N_opt,N_opt_err=curve_fit(best_line,en[[0,2]],f[[0,2]],
                              p0=b,sigma=ferr[[0,2]],absolute_sigma=True)
        N_opt_err=np.sqrt(N_opt_err[0])

        myline=lambda x: a*x+N_opt

        fe_flux=f[1]-myline(en[1])*frac
        fe_flux_err=np.sqrt(ferr[1]**2+N_opt_err**2)
        test_ratio.append(myline(en[1])*frac/(f[2]*frac))
        #fe_flux_err=np.sqrt(ferr[1]**2)

        diff.append(fe_flux)
        diff_err.append(fe_flux_err)

    diff=np.asarray(diff)
    diff=np.reshape(diff,N)
    diff_err=np.asarray(diff_err)
    diff_err=np.reshape(diff_err,N)
    return diff,diff_err,test_ratio

frac=1
fe_flux,fe_flux_err,ra=find_fe_flux(lc46, lc67, lc79,frac=frac)

# Flux_tot=lc67.rate
# Flux_fe=fe_flux
# Flux_cont=Flux_tot-Flux_fe

# eqw=(Flux_tot-Flux_cont)*dE[1]/Flux_cont*1000


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

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))



print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')
print(f'Mean Flux: {np.mean(fe_flux)}')
print(f'Std Flux: {np.std(fe_flux)}')
wmean,wstd=weighted_avg_and_std(fe_flux, fe_flux_err**(-2))
print(f'Weighted mean flux: {wmean}')
print(f'Weighted std flux: {wstd}')
print(f'Weighted mean/weighted std : {wmean/wstd}')

print(f'Mean Error: {np.mean(fe_flux_err)}')
print(f'Mean significance (mean flux/err): {np.mean(fe_flux/fe_flux_err)}')


#%% iron line flux save in fits

os.chdir(path_to_lc+'/lc712_1')
os.system(f'cp lc712AB_sr.lc_bary_orb_corr fe_line_{frac}.lc_bary_orb_corr')
with fits.open(f'/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_1/fe_line_{frac}.lc_bary_orb_corr', mode='update') as hdul:
    hdul[1].data['rate']=fe_flux
    hdul[1].data['error']=fe_flux_err
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


ccf=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_1/ccf_{frac}.qdp',skip_header=3)
N=int(ccf[:,0].shape[0]/2)
ax.errorbar(ccf[:,0],ccf[:,2],ccf[:,3],drawstyle='steps-mid')

#ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')
ax.errorbar(-ccf[:N+1,0],ccf[:N+1,2],ccf[:N+1,3],alpha=0.5,color='r',drawstyle='steps-mid')

#ax.set_xlabel('Iron line leads <---Delay, s--->  7-12 keV Flux leads')
ax.set_xlabel('Delay, s')

ax.set_ylabel('CCF')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()
plt.savefig(f'ccf_{frac}.png')



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
