#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:24:00 2020

@author: s.bykov
"""
from PipelineNuSTAR.core import *
#from Miscellaneous import pd_to_latex

matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95

import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')


ObsList=['80102002002','80102002004','80102002006','80102002008',
         '80102002010','90202031002','90202031004'] #all obsevations

savepath=f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/'



#%% pd dataframe

do_not_start_this_chunk

if input('start  calculation from the beginning?')=='y':
    ObsParams=pd.DataFrame()
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))

        nu_obs=NustarObservation(ObsID)
        nu_obs.scan_spe_results()
        ser=nu_obs.pandas_series()
        ObsParams=ObsParams.append(ser)



name='standard_nuproducts'
pd.to_pickle(ObsParams,f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/pandas_data/{name}.pkl')
ObsParams.to_csv(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/{name}.csv',index=0)

#%% load pd dataframe

results_path='/Users/s.bykov/work/xray_pulsars/nustar/plots_results/pandas_data/'

filename='standard_nuproducts'
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')



#%% plot bat maxi and nustar iron flux
fig, ax_bat = plt.subplots()
ax=ax_bat.twinx()

time=ObsParams.MJD_START
fl,fl_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_flux_gaussian_4_12',
                          funct=lambda x: 10**x/1e-10)

ax.errorbar(time,fl,fl_err,fmt='.',color='g',marker='s',ms=5)

ax.set_ylabel('NuSTAR Iron line Flux \n $10^{-10}$ cgs',color='g')
ax.set_xlabel('Time, MJD')


maxilc=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/data/J0334+531_00055054g_lc_1day_all.dat')
maxidays=maxilc[:,0]
maxirate=maxilc[:,5]
maxierror=maxilc[:,6]



time=ObsParams.MJD_START
ax_bat.errorbar(maxidays,maxirate,maxierror,fmt='.',color='b',alpha=0.7)
ax_bat.set_xlim(57180,57320)
ax_bat.set_ylim(-0.001,1.5)

ax_bat.set_ylabel('MAXI rate (4-10 keV) \n photons/s/cm^2',color='b')
ax_bat.set_xlabel('Time, MJD')


ax.legend()
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'bat_lc.pdf',dpi=500)

plt.show()



#%% plot bat maxi and nustar flux
fig, ax_bat = plt.subplots()
ax=ax_bat.twinx()

time=ObsParams.MJD_START
fl,fl_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_flux_gabslog_12_79',
                          funct=lambda x: 10**x/1e-9)

ax.errorbar(time,fl,fl_err,fmt='.',color='g',marker='s',ms=5)

ax.set_ylabel('NuSTAR Flux (12-79 keV), \n $10^{-9}$ cgs',color='g')
ax.set_xlabel('Time, MJD')

batlc=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/data/V0332p53.lc.txt')
batdays=batlc[:,0]
batrate=batlc[:,1]
baterror=batlc[:,2]

maxilc=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/data/J0334+531_00055054g_lc_1day_all.dat')
maxidays=maxilc[:,0]
maxirate=maxilc[:,5]
maxierror=maxilc[:,6]



time=ObsParams.MJD_START
#ax_bat.errorbar(batdays,batrate,baterror,fmt='.',color='b',alpha=0.7)
ax_bat.errorbar(maxidays,maxirate,maxierror,fmt='.',color='b',alpha=0.7)
#ax.set_yscale('log')
ax_bat.set_xlim(57180,57320)
ax_bat.set_ylim(-0.001,1.5)

#ax_bat.set_ylabel('BAT rate (15-50 keV) \n counts/s/cm^2',color='b')
ax_bat.set_ylabel('MAXI rate (4-10 keV) \n photons/s/cm^2',color='b')
ax_bat.set_xlabel('Time, MJD')


ax.legend()
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'bat_lc.pdf',dpi=500)

plt.show()


#%% plot eqw stuff

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_eqw_gaussian',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,label='Sigma=0.3 keV',alpha=0.8)



#eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_eqw_gaussian',funct=lambda x: 1e3*x)

#ax.errorbar(time,eqw,eqw_err,fmt='.',color='y',marker='s',ms=4,label='Sigma free',alpha=0.8)

#ax.legend(loc='best')
ax.set_xlim(57200,57320)
ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')


from Misc.doppler_correction import orb_params_v0332
T_p=orb_params_v0332['T_p']/86400
P=orb_params_v0332['P']/86400


for i in range(1,15):
    ax.axvline(T_p+i*P,color='k',ls=':',alpha=0.3)



fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw.pdf',dpi=500)

plt.show()


#%%plot iron line stuff

for par in ['Sigma3','LineE2']:
    fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

    time=ObsParams.MJD_START

    eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_'+par,funct=lambda x: x)

    ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,label='Sigma=0.3 keV',alpha=0.8)



    eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_'+par,funct=lambda x: x)

    ax.errorbar(time,eqw,eqw_err,fmt='.',color='y',marker='s',ms=4,label='Sigma free',alpha=0.8)

    ax.legend(loc='best')
    ax.set_xlim(57200,57320)
    ax.set_ylabel('Iron line '+par ,color='k')
    ax.set_xlabel('Time, MJD')
    plt.show()

#%% ============LATEX TABLES ============
from Misc.TeX_Tables import pandas_to_tex
from Misc.TeX_Tables.pandas_to_tex import *


def tex_line():
    null=lambda x: x
    free_columns=['ObsID',model+'MJD',model+'exposure',model+'chi2_red']
    free_columns_functions=[null,null,lambda x: x/1000,null]
    free_columns_formats=[0,1,0,2]

    err_columns=['LineE2','Sigma3','norm4','eqw_gaussian',
                 ]
    err_functions=[null,null,lambda x: 1000*x, lambda x: 1000*x,
                   ]
    err_formats=[2,2,1,0]

    err_columns=[model+item for item in err_columns]

    headers=['ObsID','Time, MJD','Exposure, ks','$\chi^2_{red}$',
                  'Iron Line: energy','Iron Line: sigma','Iron Line: norm','Iron Line: eq. width',
                 ]



    transpose=1
    df_tex=make_latex_table(df,
                          free_columns=free_columns, free_columns_functions=free_columns_functions,
                          free_columns_formats=free_columns_formats,
                          err_columns=err_columns, err_functions=err_functions,
                          err_formats=err_formats)
    df_tex.columns=headers
    save_latex_table(df_tex, savepath=savepath+'/tex/linepars_'+name,
                     columns_to_write='DEFAULT',
                     columns_names='DEFAULT',
                     transpose=transpose)




def tex_all():
    null=lambda x: x
    free_columns=['ObsID',model+'MJD',model+'exposure',model+'chi2_red']
    free_columns_functions=[null,null,lambda x: x/1000,null]
    free_columns_formats=[0,1,0,2]

    err_columns=['D5','Ecycle6','sigma7',
                 'D8','Ecycle9','sigma10',
                 'T012','kT13','taup14',
                 'factor17']
    err_functions=[null,null,null,
                   null,null,null,
                   null,null,null,
                   null]
    err_formats=[3,2,1,
                 3,2,1,
                 2,2,1,
                 3]

    err_columns=[model+item for item in err_columns]

    headers=['ObsID','Time, MJD','Exposure, ks','$\chi^2_{red}$',
                 'CRSF: depth', 'CRSF: energy','CRSF : sigma',
                 'CRSF (1h): depth', 'CRSF (1h): energy','CRSF (1h) : sigma',
                 'compTT: $T_0$', 'compTT: kT', 'compTT: $\tau$',
                 'C(FPMA/FPMB)']



    transpose=1
    df_tex=make_latex_table(df,
                          free_columns=free_columns, free_columns_functions=free_columns_functions,
                          free_columns_formats=free_columns_formats,
                          err_columns=err_columns, err_functions=err_functions,
                          err_formats=err_formats)
    df_tex.columns=headers
    save_latex_table(df_tex, savepath=savepath+'/tex/allpars_'+name,
                     columns_to_write='DEFAULT',
                     columns_names='DEFAULT',
                     transpose=transpose)



#%% latex tables, sigma free

indexNames = pd.core.indexes.base.Index(['obs80102002008', 'obs90202031002', 'obs90202031004'])
df=ObsParams.drop(indexNames , inplace=False)
model='comptt_gabslog_'
ParList=[col for col in df.columns if 'bin' not in col and model in col and 'sigma' not in col and '_hi' not in col and '_lo' not in col  ]
name='free_sigma.tex'
tex_line()
tex_all()


#%% latex tables, sigma fix

indexNames = pd.core.indexes.base.Index(['obs90202031002', 'obs90202031004'])
df=ObsParams.drop(indexNames , inplace=False)

model='comptt_gabslog_sigma03_'
ParList=[col for col in df.columns if 'bin' not in col and model in col and '_hi' not in col and '_lo' not in col  ]
name='fix_sigma.tex'
tex_line()
tex_all()





#%% crosscor 802

data=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc67/ccf67vs58.dat',
                   skip_header=3)
lag,ccf,err=data[:,0],data[:,2],data[:,3]
plt.errorbar(lag,ccf,err)

plt.ylim(-0.02,0.07)
plt.xlim(-3,3)



N=len(ccf)

ccf_right=ccf[N//2+1:]
lag_left=lag[:N//2]

plt.plot(lag_left,ccf_right[::-1],color='r')
plt.show()




#%% power spectrum fun

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6/2))
ax.set_xlabel('Frequency, Hz',fontsize=8)
ax.set_ylabel('Power, (RMS/mean)^2/Hz')

ObsList=['80102002002',
 '80102002004',
 '80102002006',
 '80102002008',
 '80102002010']

for ID in ['80102002002','80102002004','80102002006']:


    powspec_folder=f'/Users/s.bykov/work/xray_pulsars/nustar/results/out{ID}/products/spe_and_lc/powspec/'
    powspec_data=np.genfromtxt(powspec_folder+'powspec_AB.qdp',skip_header=3,usecols=(0,1,2,3))

    f,df,P,dP=powspec_data[:,0],powspec_data[:,1],powspec_data[:,2],powspec_data[:,3]

    P=P-np.mean(P[-5:-1])

    ax.loglog(f,P,'.',label=ID)
    ax.errorbar(f,P,dP,df,fmt='none',color='gray',alpha=0.5)

plt.legend()
plt.show()
plt.savefig(savepath+f'pdf_first.pdf',dpi=500)





fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6/2))
ax.set_xlabel('Frequency, Hz',fontsize=8)
ax.set_ylabel('Power, (RMS/mean)^2/Hz')


ObsList=['80102002002',
 '80102002004',
 '80102002006',
 '80102002008',
 '80102002010']

for ID in ['80102002008','80102002010']:


    powspec_folder=f'/Users/s.bykov/work/xray_pulsars/nustar/results/out{ID}/products/spe_and_lc/powspec/'
    powspec_data=np.genfromtxt(powspec_folder+'powspec_AB.qdp',skip_header=3,usecols=(0,1,2,3))

    f,df,P,dP=powspec_data[:,0],powspec_data[:,1],powspec_data[:,2],powspec_data[:,3]

    P=P-np.mean(P[-5:-1])

    ax.loglog(f,P,'.',label=ID)
    ax.errorbar(f,P,dP,df,fmt='none',color='gray',alpha=0.5)

plt.legend()
plt.show()
plt.savefig(savepath+f'pdf_second.pdf',dpi=500)





#%% power spectrum fun - 7-12 keV

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6/2))
ax.set_xlabel('Frequency, Hz',fontsize=8)
ax.set_ylabel('Power, (RMS/mean)^2/Hz')

ObsList=['80102002002',
 '80102002004',
 '80102002006',
 '80102002008',
 '80102002010']

for ID in ['80102002002','80102002004']:


    powspec_folder=f'/Users/s.bykov/work/xray_pulsars/nustar/results/out{ID}/products/lc712_0.1/powspec/'
    powspec_data=np.genfromtxt(powspec_folder+'powspecAB.qdp',skip_header=3,usecols=(0,1,2,3))

    f,df,P,dP=powspec_data[:,0],powspec_data[:,1],powspec_data[:,2],powspec_data[:,3]

    P=P-np.mean(P[-5:-1])

    ax.loglog(f,P,'-',label=ID,drawstyle='steps-mid')
    ax.errorbar(f,P,dP,df,fmt='none',color='gray',alpha=0.5)

plt.legend()
plt.xlim(0.0005,1)
plt.ylim(1e-3,6)
plt.show()
plt.savefig(savepath+f'pdf_first.pdf',dpi=500)





fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6/2))
ax.set_xlabel('Frequency, Hz',fontsize=8)
ax.set_ylabel('Power, (RMS/mean)^2/Hz')


ObsList=['80102002002',
 '80102002004',
 '80102002006',
 '80102002008',
 '80102002010']

for ID in ['80102002006','80102002008','80102002010']:


    powspec_folder=f'/Users/s.bykov/work/xray_pulsars/nustar/results/out{ID}/products/spe_and_lc/powspec/'
    powspec_data=np.genfromtxt(powspec_folder+'powspec_AB.qdp',skip_header=3,usecols=(0,1,2,3))

    f,df,P,dP=powspec_data[:,0],powspec_data[:,1],powspec_data[:,2],powspec_data[:,3]

    P=P-np.mean(P[-5:-1])

    ax.loglog(f,P,'-',label=ID,drawstyle='steps-mid')
    ax.errorbar(f,P,dP,df,fmt='none',color='gray',alpha=0.5)

plt.legend()
plt.xlim(0.0005,1)
plt.ylim(1e-3,6)
plt.show()
plt.savefig(savepath+f'pdf_second.pdf',dpi=500)
