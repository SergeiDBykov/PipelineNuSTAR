#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:24:00 2020

@author: s.bykov
"""

#%% imports and definitions
from PipelineNuSTAR.core import *
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')


ObsList=['80102002002','80102002004','80102002006','80102002008',
         '80102002010','90202031002','90202031004'] #all obsevations


#%% select ObsID
ObsID=ObsList[4]
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

# nu_obs.fit_spe(spe_folder='spe_and_lc',result_name='comptt_gabslog_sigma02',model='comptt_gabslog',lowlim='4',uplim='79',
#                xspec_comms=['newpar 3 0.2 -1'])

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

#%% lc for iron line timeseries
#eergy - channel
#4 kev - 60
#6 kev - 110
#7 kev - 135
#9 kev - 185
#12 keV - 260

#4-6 keV
binsize=1
for mode in ['A','B']:
    nu_obs.make_lc(mode=mode,outdir='lc46_'+str(binsize),stemout='lc46'+mode,pilow='60',pihigh='110',binsize=binsize)
#lcmath
#6-7 keV
for mode in ['A','B']:
    nu_obs.make_lc(mode=mode,outdir='lc67_'+str(binsize),stemout='lc67'+mode,pilow='110',pihigh='135',binsize=binsize)

#7-9 keV
for mode in ['A','B']:
    nu_obs.make_lc(mode=mode,outdir='lc79_'+str(binsize),stemout='lc79'+mode,pilow='135',pihigh='185',binsize=binsize)

#7-12
for mode in ['A','B']:
    nu_obs.make_lc(mode=mode,outdir='lc712_'+str(binsize),stemout='lc712'+mode,pilow='135',pihigh='260',binsize=binsize)


#orb corrs

nu_obs.orb_correction_lc(folder='lc46_'+str(binsize),filename=f'lc46AB_sr.lc_bary',q=0)
nu_obs.orb_correction_lc(folder='lc67_'+str(binsize),filename=f'lc67AB_sr.lc_bary',q=0)
nu_obs.orb_correction_lc(folder='lc79_'+str(binsize),filename=f'lc79AB_sr.lc_bary',q=0)
nu_obs.orb_correction_lc(folder='lc712_'+str(binsize),filename=f'lc712AB_sr.lc_bary',q=0)


#%% lc in 4-12 keV  range

for mode in ['A']:
    nu_obs.make_lc(mode=mode,outdir='lc412',stemout='lc412'+mode,pilow='60',pihigh='260')

nu_obs.orb_correction_lc(folder='lc412',filename='lc412A_sr.lc_bary')



#%% cross corr eqw  stuff:
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

if ObsID in ['80102002002','80102002004','80102002006']:
    #model='comptt_gabslog_sigma_free_temp_fix'
    #fluxpar='flux_gabslog_7_12'

    matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=2)
    plt.subplots_adjust(hspace=1)

    model='cutoffpl_sigma_03'
    fluxpar='flux_cutoffpl_4_12'
    nu_obs=NustarObservation(ObsID)
    nu_obs.scan_spe_results()
    ser=nu_obs.pandas_series()
    fig = plt.figure()
    rows=7
    cols=3
    ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
    ax_flux=ax_eqw.twinx()

    ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    ax_flux2=ax_fe_norm.twinx()


    ax_ccf = plt.subplot2grid((rows,cols), (5, 0), rowspan=2, colspan=3)


    #ax_eqw.set_title(nu_obs.ObsID+f'\n  model: {model}')
    phase,eqw,_=nu_obs.ph_res_param(model=model,param='eqw_gauss',funct=lambda x: 1000*x,
                      ax=ax_eqw,color='c',alpha=0.6)
    phase,norm,_=nu_obs.ph_res_param(model=model,param='norm4',funct=lambda x: x,
                      ax=ax_fe_norm,color='r',alpha=0.6)

    phase,flux712,_=nu_obs.ph_res_param(model=model,param=fluxpar,funct=lambda x: 10**x/1e-8,
                              ax=ax_flux,color='k',ls=':',alpha=0.6)
    phase,flux712,_=nu_obs.ph_res_param(model=model,param=fluxpar,funct=lambda x: 10**x/1e-8,
                              ax=ax_flux2,color='k',ls=':',alpha=0.6)

    ax_eqw.set_ylabel('Iron line \n Eq. width, eV',color='c',fontsize=8)
    ax_fe_norm.set_ylabel('Iron line \n Norm. ',color='r',fontsize=8)

    ax_flux.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_flux2.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)

    ax_fe_norm.set_xlabel('Phase',fontsize=8)
    ax_eqw.set_title(nu_obs.ObsID)

    period=ser['period_spe_and_lcA_sr.bary_lc_orb_corr']
    CCF=cross_correlation.CrossCorrelation(phase*period,eqw,flux712,circular=True)
    lag,ccf=CCF.calc_ccf()
    peaks,_,_=CCF.find_max()
    delay=min(peaks[peaks>0])
    #self.write_to_obs_info(self.fasebin_info_file,'deltat',delay)
    #self.write_to_obs_info(self.fasebin_info_file,'deltat_err',period/nph)
    #ax_ccf.axvline(delay,ls=':',color='g',alpha=0.5)

    ax_ccf.plot(lag,ccf,color='b',alpha=0.6)
    #ax_ccf.set_title(f'Flux lags <--- 0 ---> Eqw lags',fontsize=8)
    ax_ccf.set_xlim(-period,+period)
    ax_ccf.set_xlabel('Eqw Delay, sec')
    ax_ccf.set_ylabel('Pearson r')

    plt.show()

    fig.tight_layout()
    sns.despine(fig,top=1,right=0)
    savepath=f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/'
    plt.savefig(savepath+f'ph_res_{nu_obs.ObsID}_{model}_report.pdf',dpi=500)



#%% cross corr eqw  stuff 712 flux:
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

if ObsID in ['80102002002','80102002004','80102002006']:
    #model='comptt_gabslog_sigma_free_temp_fix'
    #fluxpar='flux_gabslog_7_12'

    matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=2)
    plt.subplots_adjust(hspace=1)

    model='cutoffpl_sigma_03'
    fluxpar='flux_cutoffpl_7_12'
    nu_obs=NustarObservation(ObsID)
    nu_obs.scan_spe_results()
    ser=nu_obs.pandas_series()
    fig = plt.figure()
    rows=7
    cols=3
    ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
    ax_flux=ax_eqw.twinx()

    ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    ax_flux2=ax_fe_norm.twinx()


    ax_ccf = plt.subplot2grid((rows,cols), (5, 0), rowspan=2, colspan=3)


    #ax_eqw.set_title(nu_obs.ObsID+f'\n  model: {model}')
    phase,eqw,_=nu_obs.ph_res_param(model=model,param='eqw_gauss',funct=lambda x: 1000*x,
                      ax=ax_eqw,color='c',alpha=0.6)
    phase,norm,_=nu_obs.ph_res_param(model=model,param='norm4',funct=lambda x: x,
                      ax=ax_fe_norm,color='r',alpha=0.6)

    phase,flux712,_=nu_obs.ph_res_param(model=model,param=fluxpar,funct=lambda x: 10**x/1e-8,
                              ax=ax_flux,color='k',ls=':',alpha=0.6)
    phase,flux712,_=nu_obs.ph_res_param(model=model,param=fluxpar,funct=lambda x: 10**x/1e-8,
                              ax=ax_flux2,color='k',ls=':',alpha=0.6)

    ax_eqw.set_ylabel('Iron line \n Eq. width, eV',color='c',fontsize=8)
    ax_fe_norm.set_ylabel('Iron line \n Norm. ',color='r',fontsize=8)

    ax_flux.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_flux2.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)

    ax_fe_norm.set_xlabel('Phase',fontsize=8)
    ax_eqw.set_title(nu_obs.ObsID)

    period=ser['period_spe_and_lcA_sr.bary_lc_orb_corr']
    CCF=cross_correlation.CrossCorrelation(phase*period,eqw,flux712,circular=True)
    lag,ccf=CCF.calc_ccf()
    peaks,_,_=CCF.find_max()
    delay=min(peaks[peaks>0])
    #self.write_to_obs_info(self.fasebin_info_file,'deltat',delay)
    #self.write_to_obs_info(self.fasebin_info_file,'deltat_err',period/nph)
    #ax_ccf.axvline(delay,ls=':',color='g',alpha=0.5)

    ax_ccf.plot(lag,ccf,color='b',alpha=0.6)
    #ax_ccf.set_title(f'Flux lags <--- 0 ---> Eqw lags',fontsize=8)
    ax_ccf.set_xlim(-period,+period)
    ax_ccf.set_xlabel('Eqw Delay, sec')
    ax_ccf.set_ylabel('Pearson r')

    plt.show()

    fig.tight_layout()
    sns.despine(fig,top=1,right=0)
    savepath=f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/'
    plt.savefig(savepath+f'ph_res_{nu_obs.ObsID}_{model}_report.pdf',dpi=500)




#%% lc in 7-12 keV range for the power spectra


#eergy - channel
#4 kev - 60
#6 kev - 110
#7 kev - 135
#9 kev - 185
#12 keV - 260

binsize=0.1

#7-12
for mode in ['A','B']:
    nu_obs.make_lc(mode=mode,outdir='lc712_'+str(binsize),stemout='lc712'+mode,pilow='135',pihigh='260',binsize=binsize)


#orb corrs

# nu_obs.orb_correction_lc(folder='lc46_'+str(binsize),filename=f'lc46AB_sr.lc_bary',q=0)
# nu_obs.orb_correction_lc(folder='lc67_'+str(binsize),filename=f'lc67AB_sr.lc_bary',q=0)
# nu_obs.orb_correction_lc(folder='lc79_'+str(binsize),filename=f'lc79AB_sr.lc_bary',q=0)
#nu_obs.orb_correction_lc(folder='lc712_'+str(binsize),filename=f'lc712AB_sr.lc_bary',q=0)
