#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:24:00 2020

@author: s.bykov
"""
from PipelineNuSTAR.core import *
from Miscellaneous import pd_to_latex
import seaborn as sns
sns.set_palette("pastel")

ObsList=['80102002002','80102002004','80102002006','80102002008',
         '80102002010','90202031002','90202031004'] #all obsevations




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



#%% plot bat flux


batlc=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/data/V0332p53.lc.txt')
savepath=f'/Users/s.bykov/work/xray_pulsars/nustar/plots_results/'
batdays=batlc[:,0]
batrate=batlc[:,1]
baterror=batlc[:,2]

matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.right']=0.95

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})
# ax_bat=ax.twinx()
# ax_bat,ax=ax,ax_bat
# time=ObsParams.MJD_START
# fl,fl_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma02_flux_gabslog_12_79',funct=lambda x: 10**x/1e-9)
ax_bat.errorbar(time,fl,fl_err,fmt='.',color='g')
ax.errorbar(batdays,batrate,baterror,fmt='.',alpha=0.7)
#ax.set_yscale('log')
ax.set_xlim(57180,57610)
ax.set_ylim(-0.001,0.25)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('BAT rate (15-50 keV) \n counts/s/cm^2')

fig.savefig(savepath+f'bat_lc.png')
plt.show()


#%% plot eqw stuff



eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma02_eqw_gaussian',funct=lambda x: 1e3*x)

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})
ax_bat=ax.twinx()
ax.errorbar(time,eqw,eqw_err,fmt='.',label='eqw_sigma_0.2',alpha=0.7)

ax_bat.errorbar(batdays,batrate,baterror,fmt='.',color='g',alpha=0.7)
ax.set_xlim(57180,57610)


plt.show()



time=ObsParams.MJD_START
eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_eqw_gaussian',funct=lambda x: 1e3*x)

#fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})
ax_bat=ax.twinx()
ax.errorbar(time,eqw,eqw_err,fmt='.',color='m',label='eqw_sigma_free',alpha=0.7)

ax.set_ylabel('Eq width, eV')
ax_bat.set_ylim(-0.001,0.25)
ax_bat.set_xlabel('Time, MJD')
ax_bat.set_ylabel('BAT rate (15-50 keV) \n counts/s/cm^2')

from Miscellaneous.doppler_correction import orb_params_v0332
T_p=orb_params_v0332['T_p']/86400
P=orb_params_v0332['P']/86400


for i in range(1,15):
    ax.axvline(T_p+i*P,color='k',ls=':',alpha=0.3)


ax.legend()
fig.savefig(savepath+f'eqw.png')
plt.show()





#%% latex table constant sigma
from Miscellaneous import pd_to_latex as latex
transpose=1
model='comptt_gabslog_sigma02_'
la_pd=latex.make_latex_table(ObsParams,
                 ['MJD_START','ObsID',model+'chi2_red'],
                 [lambda x:x,lambda x:x,lambda x:x],
                 [1,1,2],
                 [model+'eqw_gaussian',model+'flux_gabslog_4_79'],
                 [lambda x: 1000*x,lambda x: 10**(x)/1e-9],
                 [0,3],
                 transpose=transpose)


la_pd.to_latex(escape = False, index = transpose)
la_pd.to_latex(buf=savepath+'table_sigma02.tex',escape = False, index = 1)


from Miscellaneous import pd_to_latex as latex
transpose=1
model='comptt_gabslog_'
la_pd=latex.make_latex_table(ObsParams,
                 ['MJD_START','ObsID',model+'chi2_red'],
                 [lambda x:x,lambda x:x,lambda x:x],
                 [1,1,2],
                 [model+'eqw_gaussian',model+'Sigma3'],
                 [lambda x: 1000*x,lambda x: x],
                 [0,3],
                 transpose=transpose)


la_pd.to_latex(escape = False, index = transpose)
la_pd.to_latex(buf=savepath+'table_sigma_free.tex',escape = False, index = 1)



#%% latex ultra table
nof=lambda x:x
model='comptt_gabslog_' #'comptt_gabslog_sigma02_' or 'comptt_gabslog_'
# free_columns=['exposure','chi2_red','LineE2','Sigma3','norm4','eqw_gaussian',
#               'D5','Ecycle6','sigma7','D8','Ecycle9','sigma10',
#               'T012','kT13','taup14','norm16','factor17'
#               ]
# free_columns_functions=[nof,nof,nof,nof,lambda x: 1000*x,lambda x: 1000*x,
#               nof,nof,nof,nof,nof,nof,
#               nof,nof,nof,nof,nof]
# free_columns_formats=[0,3,2,2,2,1,1,
#               3,2,1,3,2,1,
#               2,2,1,3,3]


free_columns=['exposure','chi2_red']
free_columns_functions=[nof,nof]
free_columns_formats=[0,3]

free_columns=[model+item for item in free_columns]


err_columns=['LineE2','Sigma3','norm4','eqw_gaussian',
              'D5','Ecycle6','sigma7','D8','Ecycle9','sigma10',
              'T012','kT13','taup14','norm16','factor17'
              ]
err_functions=[nof,nof,lambda x: 1000*x,lambda x: 1000*x,
              nof,nof,nof,nof,nof,nof,
              nof,nof,nof,lambda x: x*1000,nof]
err_formats=[2,2,2,0,
             3,3,2,3,3,2,
             2,2,1,0,3]
err_columns=[model+item for item in err_columns]



transpose=1
la_pd=latex.make_latex_table(ObsParams,
                             free_columns=free_columns, free_columns_functions=free_columns_functions, free_columns_formats=free_columns_formats,
                             err_columns=err_columns, err_functions=err_functions, err_formats=err_formats,
                             transpose=transpose)

tablename=model
la_pd.to_latex(escape = False, index = transpose)
la_pd.to_latex(buf=savepath+tablename+'.tex',escape = False, index = transpose)
