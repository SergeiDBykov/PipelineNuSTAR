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



#%% plot bat and nustar flux
fig, ax_bat = plt.subplots()
ax=ax_bat.twinx()

time=ObsParams.MJD_START
fl,fl_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_flux_gabslog_12_79',
                          funct=lambda x: 10**x/1e-9)

ax.errorbar(time,fl,fl_err,fmt='.',color='g',marker='s',ms=5)

ax.set_ylabel('NuSTAR Flux (12-79 keV), \n 10^(-9) cgs',color='g')
ax.set_xlabel('Time, MJD')

batlc=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/nustar/data/V0332p53.lc.txt')
batdays=batlc[:,0]
batrate=batlc[:,1]
baterror=batlc[:,2]



time=ObsParams.MJD_START
ax_bat.errorbar(batdays,batrate,baterror,fmt='.',color='b',alpha=0.7)
#ax.set_yscale('log')
ax_bat.set_xlim(57180,57320)
ax_bat.set_ylim(-0.001,0.25)

ax_bat.set_ylabel('BAT rate (15-50 keV) \n counts/s/cm^2',color='b')

ax.legend()
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'bat_lc.png',dpi=500)

plt.show()


#%% plot eqw stuff

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_sigma03_eqw_gaussian',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,label='Sigma=0.3 keV',alpha=0.8)



eqw,eqw_err=vals_and_errors(ObsParams,'comptt_gabslog_eqw_gaussian',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='y',marker='s',ms=4,label='Sigma free',alpha=0.8)

ax.legend(loc='best')
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
plt.savefig(savepath+f'eqw.png',dpi=500)

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

#indexNames = ObsParams[ ObsParams['ObsID'] > 80102002007 ].index
#df=ObsParams.drop(indexNames , inplace=False)
df=ObsParams
model='comptt_gabslog_sigma03_'
ParList=[col for col in df.columns if 'bin' not in col and model in col and '_hi' not in col and '_lo' not in col  ]
name='fix_sigma.tex'
tex_line()
tex_all()


#%% trash
'''
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
'''