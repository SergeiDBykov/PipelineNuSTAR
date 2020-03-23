#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:24:00 2020

@author: s.bykov
"""

#%% imports and definitions
from PipelineNuSTAR.core import *

ObsList=['80102002002','80102002004','80102002006','80102002008',
         '80102002010','90202031002','90202031004'] #all obsevations


#%% select ObsID
ObsID=ObsList[2]
nu_obs=NustarObservation(ObsID)

STOP
#%% pipeline and region  check
nu_obs.nupipeline()
nu_obs.make_regions()

#%% basic nuproducts
for mode in ['A','B']:
    nu_obs.nuproducts(outdir='spe_and_lc',stemout='spe_and_lc'+mode,mode=mode)

#%% spectral fit with complex model

nu_obs.fit_spe(spe_folder='spe_and_lc',result_name='comptt_gabslog',model='comptt_gabslog',lowlim='4',uplim='79')

nu_obs.fit_spe(spe_folder='spe_and_lc',result_name='comptt_gabslog_sigma02',model='comptt_gabslog',lowlim='4',uplim='79',
               xspec_comms=['newpar 3 0.2 -1'])

nu_obs.fit_spe(spe_folder='spe_and_lc',result_name='comptt_gabslog_sigma03',model='comptt_gabslog',lowlim='4',uplim='79',
               xspec_comms=['newpar 3 0.3 -1'])

#%% write fit results to spe info file
nu_obs.scan_spe_results()


#%% period finding: orb correction
for mode in ['A','B']:
    nu_obs.orb_correction_lc(filename=f'spe_and_lc{mode}_sr.bary_lc',q=0)


#%% period finding: efsearch
nu_obs.make_efsearch('spe_and_lc')


#%% period finding: efsearch fit
nu_obs.find_period('spe_and_lc')




# #%% make gtis and create ph_res_script
# for mode in ['A','B']:
#     ser=nu_obs.pandas_series()
#     per_col=[col for col in ser.index.values if 'period' in col and 'lc'+mode in col and 'err' not in col]
#     period_mode=ser[per_col].values[0]
#     period_mode='%.4f'% period_mode
#     print('PERIOD '+mode, period_mode)
# period_mode=input('Enter period value for phase_resolved spectroscopy: ')
# for mode in ['A','B']:
#         nu_obs.make_gti_from_lc('spe_and_lc',mode=mode,period=period_mode)
#         nu_obs.phase_resolved_spectra(mode=mode)

#%% fasebin for observation
nu_obs.fasebin()

#%% ph_res_ orb_ corr of lc
os.chdir(nu_obs.products_path)
os.chdir('phase_resolved')
files=glob(f'*_sr*bary*')
for filename in files:
    nu_obs.orb_correction_lc(folder='phase_resolved',filename=filename)


#%% ph_res fit
nu_obs.scan_spe_results()
ser=nu_obs.pandas_series()
kT=ser.comptt_gabslog_kT13
T0=ser.comptt_gabslog_T012

# xspec_comm_g=[f'newpar 12 {T0} -1', f'newpar 13 {kT} -1','newpar 3 0.2 -1']
# nu_obs.fit_ph_res_spe(result_name='comptt_gabslog_sigma02', model='comptt_gabslog',
#                       lowlim=4,uplim=79,
#                       xspec_comms=xspec_comm_g)

xspec_comm_g=[f'newpar 12 {T0} -1', f'newpar 13 {kT} -1','newpar 3 0.3']
nu_obs.fit_ph_res_spe(result_name='comptt_gabslog_sigma_free_temp_fix', model='comptt_gabslog',
                      lowlim=4,uplim=79,
                      xspec_comms=xspec_comm_g)


xspec_comm_g=['newpar 3 0.3']
nu_obs.fit_ph_res_spe(result_name='comptt_gabslog_sigma_free_temp_free', model='comptt_gabslog',
                      lowlim=4,uplim=79,
                      xspec_comms=xspec_comm_g)


xspec_comm_g=['newpar 3 0.3 -1']
nu_obs.fit_ph_res_spe(result_name='comptt_gabslog_sigma_03_temp_free', model='comptt_gabslog',
                      lowlim=4,uplim=79,
                      xspec_comms=xspec_comm_g)


nu_obs.fit_ph_res_spe(model='cutoffpl',result_name='cutoffpl_sigma_free',lowlim=4,uplim=12,
                      xspec_comms=['newpar 3 0.3 '])

nu_obs.fit_ph_res_spe(model='cutoffpl',result_name='cutoffpl_sigma_03',lowlim=4,uplim=12,
                      xspec_comms=['newpar 3 0.3 -1'])

# xspec_comm_no_g=[f'newpar 9 {T0} -1', f'newpar 10 {kT} -1']
# nu_obs.fit_ph_res_spe(model='comptt_gabslog_no_gauss',result_name='comptt_gabslog_no_gauss',
#                       lowlim=4,uplim=79,
#                       xspec_comms=xspec_comm_no_g)


#%% ph res spectra results
nu_obs.scan_spe_results()
ser=nu_obs.pandas_series()
for model in ['comptt_gabslog_sigma_free_temp_free',
              'comptt_gabslog_sigma_free_temp_fix',
              'comptt_gabslog_sigma_03_temp_free',
              'cutoffpl_sigma_free',
              'cutoffpl_sigma_03']:
    try:
        nu_obs.ph_res_results(model=model)
        plt.show()
    except:
        print(f'model {model} not found')

#%% check efold and flux

fig,ax=plt.subplots()

#ax_efold=ax.twinx()
ax_efold=ax
try:
    phase,flux,err=nu_obs.ph_res_param(model='comptt_gabslog_sigma_free_temp_free',
                                   param='flux_gabslog_4_12',ax=ax,
                                   funct=lambda x: 10**x/np.mean(10**x),color='k',lw=2)
except:
    phase,flux,err=nu_obs.ph_res_param(model='cutoffpl_sigma_free',
                                   param='flux_cutoffpl_4_12',ax=ax,
                                   funct=lambda x: 10**x/np.mean(10**x),color='k',lw=2)

for file in [nu_obs.products_path+'/lc412/efold.dat',nu_obs.products_path+'/lc412/efold_25.dat']:
    efold_file=np.genfromtxt(file)

    efph=efold_file[:,0]
    efr=efold_file[:,2]
    eferr=efold_file[:,3]
    ax_efold.errorbar(efph,efr,yerr=eferr,drawstyle='steps-mid',color=numpy.random.rand(3,))

# # =============================================================================
# # efolds for obs 802 were displaced by halfbin: 802: 17223.41217923250(indef)+4.3764/8/2/86400 as a null phase value
# #   804  17275.94671724814+4.3759/8/2/86400-2*4.3759/8/86400 (2 bins shift)
# # =============================================================================
plt.show()

#%% lc in 4-12 keV  range

for mode in ['A']:
    nu_obs.make_lc(mode=mode,outdir='lc412',stemout='lc412'+mode,pilow='60',pihigh='260')

nu_obs.orb_correction_lc(folder='lc412',filename='lc412A_sr.lc_bary')


#%% cross corr stuff:
#%% lc in 6-7 keV  range

for mode in ['A']:
    nu_obs.make_lc(mode=mode,outdir='lc67',stemout='lc67'+mode,
                   pilow='110',pihigh='135',binsize=0.1)

#%% lc in 6-7 keV  range

for mode in ['A']:
    nu_obs.make_lc(mode=mode,outdir='lc58',stemout='lc58'+mode,
                   pilow='85',pihigh='160',binsize=0.1)


nu_obs.orb_correction_lc(folder='lc58',filename='lc58A_sr.lc_bary')



#%%
# TRASH
'''
#%%scan


nu_obs.scan_spe_results()
ser=nu_obs.pandas_series()
ph=np.arange(0,8)/8

eqw=[1000*ser[f'cutoffpl_bin{i}_eqw_gaussian'] for i in range(1,9)]
fl=[10**ser[f'cutoffpl_bin{i}_flux_cutoffpl_4_12'] for i in range(1,9)]


fig,ax=plt.subplots()

ax2=ax.twinx()

#ax.plot(ph,eqw,'r:')
ax2.plot(ph,fl,'k:o')


# efold=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc412/lc412A_sr.lc_bary_orb_corr.efold8',
#                     skip_header=3)
# eph=efold[:,0]
# er=efold[:,2]

# ax3=ax2.twinx()

# ax3.plot(eph,er,'m:o')
# ax3.plot(eph-2/8,er,'g-.o')
plt.show()

#%% zero times
#os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc412')

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out90202031004/products/phase_resolved')
period=4.3764

ff_orb=fits.open('phase_resolved_A_bin1_sr.lc_bary_orb_corr')
# mjdrefi=ff_bary[1].header['mjdrefi']
# mjdreff=ff_bary[1].header['mjdreff']
# timezero=ff_bary[1].header['timezero']
#time_orig_mjd=(ff_bary[1].data['time']+timezero)/86400+mjdreff+mjdrefi
time_orb_corr=ff_orb[1].data['time']


# ff_bary=fits.open('phase_resolved_A_bin1_sr.lc_bary')
# mjdrefi=ff_bary[1].header['mjdrefi']
# mjdreff=ff_bary[1].header['mjdreff']
# timezero=ff_bary[1].header['timezero']
# time_bary=(ff_bary[1].data['time']+timezero)+(mjdreff+mjdrefi)*86400
# diff_time=time_orb_corr[0]-time_bary[0]


ff_orig=fits.open('phase_resolved_A_bin1_sr.lc')
mjdrefi=ff_orig[1].header['mjdrefi']
mjdreff=ff_orig[1].header['mjdreff']
timezero=ff_orig[1].header['timezero']
time_orig=(ff_orig[1].data['time']+timezero)+(mjdreff+mjdrefi)*86400
diff_time=time_orb_corr[0]-time_orig[0]


ff_orb_mean_lc=fits.open('../spe_and_lc/spe_and_lcA_sr.bary_lc_orb_corr')
# mjdrefi=ff_bary[1].header['mjdrefi']
# mjdreff=ff_bary[1].header['mjdreff']
# timezero=ff_bary[1].header['timezero']
#time_orig_mjd=(ff_bary[1].data['time']+timezero)/86400+mjdreff+mjdrefi
time_orb_corr_mean_lc=ff_orb_mean_lc[1].data['time']


#%% test efolds
rate=[]
for i in range(1,9):
    name=f'phase_resolved_A_bin{i}_sr.lc_bary_orb_corr'
    tmpf=fits.open(name)
    tmp=tmpf[1].data['rate'].mean()
    rate.append(tmp)
    tmpf.close()
ph=(np.arange(1,9)-1)/8
fig,ax=plt.subplots()

ax.plot(ph,rate,'k-.')
ax.plot(ph+4/8,rate,'k:')
efoldfile='/Users/s.bykov/work/xray_pulsars/nustar/results/out90202031004/products/spe_and_lc/efoldA8.dat'

efold=np.genfromtxt(efoldfile,
                    skip_header=3)
eph=efold[:,0]
er=efold[:,2]

ax2=ax.twinx()

ax2.plot(eph,er,'g:o')

plt.show()


'''