#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:24:00 2020

@author: s.bykov
"""

#%% imports and definitions bla bla blah
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from glob import glob
import os
from Misc.TimeSeries import cross_corr
from Misc import  doppler_correction
from scipy.optimize import curve_fit
import PipelineNuSTAR.xspec  as nu_xspec
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')



from subprocess import call
#%% functions
def open_dir():
    call(['open',os.getcwd()])

def open_dir_in_term():
    import appscript
    appscript.app('Terminal').do_script(f'cd {os.getcwd()}')


plt.ioff()


Nu_path='/Users/s.bykov/work/xray_pulsars/nustar/'
xspec_scripts_path='/Users/s.bykov/work/xray_pulsars/nustar/python_nu_pipeline/xspec_scripts/'

def start_stop(a, trigger_val=1):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False,a==trigger_val,False]

    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])

    return idx.reshape(-1,2)-[0,1]

def vals_and_errors(ObsParams,name,funct=lambda x: x):
    if isinstance(ObsParams,pd.core.frame.DataFrame):
        par,Min,Max=funct(ObsParams[name].values),funct(ObsParams[name+'_lo'].values),funct(ObsParams[name+'_hi'].values)
    elif isinstance(ObsParams,pd.core.series.Series):
        par,Min,Max=funct(ObsParams[name]),funct(ObsParams[name+'_lo']),funct(ObsParams[name+'_hi'])

    low = par - Min
    hi  =  Max - par
    parr = par

    err=np.vstack((low,hi))
    #return np.array((parr,err))
    return parr,err

def sum_error(a,b,da,db):
    f=a+b
    sigma=np.sqrt(da**2+db**2)
    return f,sigma

def ratio_error(a,b,da,db):
    f=a/b
    sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
    return f, sigma



def gauss(t,t0,sigma,N):
    return N*np.exp(-(t-t0)**2/(2*sigma**2))/np.sqrt(sigma)

def find_columns_by_entry(series,keywals):
    pass



def fit_efsearch_data(efsearcf_fits_file,savefig=1,fit_data=1):
    efsearch=fits.open(efsearcf_fits_file)
    period=efsearch[1].data['period']
    chisq=efsearch[1].data['chisqrd1']
    p0=[period[np.argmax(chisq)],0.001,max(chisq)*np.sqrt(0.001)]

    fig,ax=plt.subplots()

    try:
        popt,perr=curve_fit(gauss,period,chisq,p0=p0)
        perr=np.sqrt(np.diag(perr))
        ax.plot(period,gauss(period,*p0),'r:',alpha=0.5)
        ax.plot(period,chisq)
        ax.plot(period,gauss(period,*popt),'k:')
        ax.set_title('Period='+str(popt[0])+'  sigma='+str(popt[1])+'\n')
        ax.set_xlabel('Period')
        ax.set_ylabel('\chi^2')
        if savefig:
            fig.savefig(f'{efsearcf_fits_file}.png')
            plt.close(fig)
        return popt[0],perr[0]
    except:
        ax.plot(period,chisq)
        ax.plot(period,gauss(period,*p0),'r-.')
        ax.set_title('failed fit'+'\n')
        ax.set_xlabel('Period')
        ax.set_ylabel('\chi^2')
        if savefig:
            fig.savefig(f'{efsearcf_fits_file}.png')
            plt.close(fig)
        return None,None


#def run_command(cmd,logpath,cmd_name='cmd'):
#    print('Running command:',cmd)
#    print('Writing to file: ',cmd_name)
#    os.system(f"echo '\n {cmd} | tee -a {logpath}/pipeline_log.txt' >> {cmd_name}.sh")
#    os.system(f'chmod -x {cmd_name}.sh')
#    return os.path.abspath(cmd_name)

def run_command(cmd,logpath=None,cmd_name='cmd',rewrite=0):
    print('Running command:',cmd)
    print('Writing to file: ',cmd_name)

    if rewrite:
        os.system(f'rm -f {cmd_name}.sh')

    #if logpath!=None:
    #    os.system(f"echo '\n {cmd} | tee -a {logpath}/pipeline_log.txt' >> {cmd_name}.sh")
    #else:
    if 1:
        os.system(f"echo '\n {cmd}' >> {cmd_name}.sh")
    os.system(f'chmod -x {cmd_name}.sh')
    return os.path.abspath(cmd_name)


def create_dir(dir):
    os.system(f'mkdir -p {dir}')

#%% NuSTAR functions and class
def PIfromE(E):
    return (E-1.6)/0.04

def EfromPI(PI):
    return PI*0.04+1.6

fpm='FPM'
modules=['A','B']


bright_obs=['80102002002','80102002004','80102002006','80102002008'] #when status expt is demanded

ObsList=['80102002002','80102002004','80102002006','80102002008',
         '80102002010','90202031002','90202031004'] #all obsevations

class NustarObservation():
    def __init__(self,ObsID):
        '''
        we load obs name and write pathes to self
        create if necessary path for results of analysis

        '''
        print('###')
        print(f'Observation {ObsID} loaded successfully')
        self.ObsID=ObsID
        self.data_path=Nu_path+'data/'+ObsID+'/'
        os.chdir(Nu_path+'results/')


        create_dir('out'+self.ObsID)
        os.chdir('out'+ObsID)
        out_path=os.getcwd()
        self.out_path=out_path
        create_dir('products')
        os.chdir('products')
        self.products_path=os.getcwd()
#        for folder in ['spe_and_lc',
#                       'lc_410','lc_1079','lc_712',
#                       'fasebin']:
#            create_dir(folder)


        os.chdir(self.out_path)
        create_dir('obs_info')
        os.chdir('obs_info')
        self.obs_info_file=open(f"obs_info_{ObsID}.txt","a+")
        self.obs_info_file.close()
        self.obs_info_file=os.path.abspath(f'obs_info_{ObsID}.txt')

        self.spe_info_file=open(f"spe_info_{ObsID}.txt","a+")
        self.spe_info_file.close()
        self.spe_info_file=os.path.abspath(f'spe_info_{ObsID}.txt')

        self.per_info_file=open(f"per_info_{ObsID}.txt","a+")
        self.per_info_file.close()
        self.per_info_file=os.path.abspath(f'per_info_{ObsID}.txt')

        self.fasebin_info_file=open(f"fasebin_info_{ObsID}.txt","a+")
        self.fasebin_info_file.close()
        self.fasebin_info_file=os.path.abspath(f'fasebin_info_{ObsID}.txt')

        os.chdir(self.out_path)
        self.write_to_obs_info(self.obs_info_file,'ObsID',self.ObsID)
        try:
            ff=fits.open(glob('*A01_cl.evt')[0])
            MJD=ff[1].header['mjd-obs']
            ff.close()
            self.write_to_obs_info(self.obs_info_file,'MJD_START',MJD)
        except:
            pass



    def write_to_obs_info(self,file,name,value,rewrite=0):

        '''
        store variables in file
        mode is a+ -append
        save name of variable and its value, for instance
        period 4.37
        '''
        mode='a+'
        if rewrite:
            mode='w+'
        with open(file,mode) as f:
            f.write(name+'\t'+str(value)+'\n')
            f.close()


    def pandas_series(self,read_obs_info_file=True):
        '''
        lots of those series in a dataframe give an opportunity to
        choose observation for general analysis (period finding, spectral approximation, etc)
        and for phase resolved spectroscopy, and to print datamodes for sa/se data
        '''
        obs_info={}
        if read_obs_info_file:
            for filepath in [self.obs_info_file,self.spe_info_file,
                             self.per_info_file,self.fasebin_info_file]:
                if os.stat(filepath).st_size==1: #strange bug: some files are not empty but host 1 carriage return, hence the error. i empty such files by hand.
                    print('1 byte filesize')
                    os.system(f'>{filepath}')
                with open(filepath) as f:
                    for line in f:
                        (key,val) = line.split()            #in case of normal entry as 'chi2 23' or 'se_config E123um'
                        try:
                            obs_info[key] = float(val)
                        except:
                            obs_info[key] = str(val)
                    f.close()

        Series=pd.Series(obs_info,name='obs'+self.ObsID)
        return Series


    def nupipeline(self):

        os.chdir(self.out_path)
        if len(glob('*A01_cl.evt'))!=0:
            raise Exception('nupipeline has already been launched!')



        nupipeline=f'''nupipeline indir={self.data_path} steminputs=nu{self.ObsID} outdir={self.out_path} obsmode=SCIENCE createexpomap=yes'''
        if self.ObsID in bright_obs:
            nupipeline+=' statusexpr="STATUS==b0000xxx00xxxx000"'
        run_command(nupipeline,logpath=self.out_path)

    def make_regions(self):
        ds9=f"ds9 -tile nu{self.ObsID}A01_cl.evt -scale log  -tile nu{self.ObsID}B01_cl.evt  -scale log -lock scale"

        run_command(ds9,self.out_path,'ds9_set')

        ds9=f"ds9 -tile nu{self.ObsID}A01_cl.evt -regions srcA.reg -tile nu{self.ObsID}A01_cl.evt   -regions bkgA.reg -tile nu{self.ObsID}B01_cl.evt -region srcB.reg  -tile nu{self.ObsID}B01_cl.evt -regions bkgB.reg  -scale log  -lock scale -lock frame wcs"
        run_command(ds9,self.out_path,'ds9_check')



    def nuproducts(self,outdir,mode, stemout='DEFAULT',
                   bkgextract='yes',
                   phafile='DEFAULT',bkgphafile='DEFAULT',runmkrmf='yes',runmkarf='yes',
                   lcfile='DEFAULT',bkglcfile='DEFAULT',
                   imagefile='DEFAULT',
                   pilow='60',pihigh='1935',binsize=0.01,
                   usrgtifile='NONE',usrgtibarycorr='no',
                   grppha=0,barycor=0):
        '''
        https://heasarc.gsfc.nasa.gov/lheasoft/ftools/caldb/help/nuproducts.html
        mode- A or B
        outdir- folder name to be created, like spe_and_lc, or lc1070 etc
        others parameters from nuproducts

        output parameters are needen to catch output pathes and run grppha/barycor or other functions
        on a fly

        '''
        if mode not in ['A','B']:
            raise Exception('choose module name from A,B')
        os.chdir(self.products_path)
        create_dir(outdir)


        ObsID=self.ObsID
        outdi=self.products_path+'/'+outdir

        nuproducts=f'''
    nuproducts \
    indir={self.out_path} \
    instrument={fpm}{mode} \
    steminputs=nu{ObsID} \
    stemout={stemout} \
    outdir={outdi} \
    srcregionfile={self.out_path}/src{mode}.reg \
    bkgextract={bkgextract} \
    bkgregionfile={self.out_path}/bkg{mode}.reg \
    binsize={binsize} \
    lcfile={lcfile} \
    phafile={phafile} \
    bkglcfile={bkglcfile} \
    bkgphafile={bkgphafile} \
    imagefile={imagefile} \
    usrgtifile={usrgtifile} \
    pilow={pilow} pihigh={pihigh} \
    usrgtibarycorr={usrgtibarycorr}\
    runmkarf={runmkarf}\
    runmkrmf={runmkrmf}'''
        cmd_path=run_command(nuproducts,self.out_path,outdir,rewrite=0)
        os.chdir(outdir)
        prod_path=os.getcwd()
        prod_path=os.path.abspath(prod_path)
        output=glob('*')
        print('created files:', output)
        print('cmd_path,products_path=',cmd_path,prod_path)


        out_spe_file=os.path.abspath(prod_path+f'/{stemout}_sr.pha')
        print('expected filename', out_spe_file)
        if grppha: self.grppha(out_spe_file,cmd_path=cmd_path)

        out_lc_file=os.path.abspath(prod_path+f'/{stemout}_sr.lc')
        print('expected filename', out_lc_file)
        if barycor: self.barycorr(out_lc_file, cmd_path)


        return cmd_path,prod_path

    def make_spe(self,group_min=25,**kwargs):
        '''
        version of nuproducts without LC

        '''
        cmd_path,prod_path=self.nuproducts(**kwargs,
                   lcfile='NONE',bkglcfile='NONE',imagefile='NONE')


        stemout=kwargs['stemout']
        out_spe_file=os.path.abspath(prod_path+f'/{stemout}_sr.pha')
        print('excpected filename', out_spe_file)

        self.grppha(out_spe_file,cmd_path=cmd_path,group_min=group_min)


    def grppha(self,infile,cmd_path,group_min=30):
        outfile=infile[0:-4]+'.pi'
        grppha=f'''grppha infile="{infile}" outfile="{outfile}"  comm="group min {group_min} & exit" clobber=yes'''
        run_command(grppha,self.out_path,cmd_path)




    def make_lc(self,**kwargs):
        '''
        version of nuproducts without spectra

        '''
        cmd_path,prod_path=self.nuproducts(**kwargs,
                   phafile='NONE',bkgphafile='NONE',runmkrmf='NONE',runmkarf='NONE',imagefile='NONE')

        stemout=kwargs['stemout']
        out_lc_file=os.path.abspath(prod_path+f'/{stemout}_sr.lc')
        print('excpected filename', out_lc_file)
        print(os.getcwd())
        self.barycorr(out_lc_file,cmd_path=cmd_path)


    def barycorr(self,infile,cmd_path,barytime='no'):
        outfile=infile+'_bary'
        barycorr=f'''      barycorr \
    infile={infile} \
    outfile={outfile}\
    orbitfiles={self.data_path}auxil/nu{self.ObsID}_orb.fits.gz\
    ra=53.7500 \
    dec=53.1733 \
    barytime={barytime} '''
        run_command(barycorr,self.out_path,cmd_path)




    def make_efsearch(self,lc_folder,p0='4.3763'):
        os.chdir(self.products_path+'/'+lc_folder)
        lcs=glob('*.bary_lc_orb_corr')
        if lcs==[]:
            raise Exception('No barycentred files found')
        print('Found barycentred lcs: ',lcs)

        for lcfile in lcs:
            efsearch=f'efsearch cfile1="{lcfile}" dper={p0} nphase=8 dres=0.00001 nper=128 outfile="{lcfile}.efs" window="-" sepoch=INDEF nbint=INDEF plot=yes plotdev="/xw"'
            run_command(efsearch,logpath=self.out_path,cmd_name='efsearch')


    def find_period(self,lc_folder,rewrite=1):
        '''
        efsearch files are fitted with gaussian.
        it also writes period and  period sigma to per info file

        '''
        os.chdir(self.products_path+'/'+lc_folder)
        lcs=glob('*.bary_lc_orb_corr')
        if lcs==[]:
            raise Exception('No barycentred files found')
        print('Found barycentred lcs: ',lcs)


        if rewrite:
            fname=self.per_info_file
            os.system(f'> {fname}') #remove content from per_info

        for lcfile in lcs:
            per,err=fit_efsearch_data(lcfile+'.efs')
            self.write_to_obs_info(self.per_info_file,f'period_{lcfile}',per)
            self.write_to_obs_info(self.per_info_file,f'period_err_{lcfile}',err)


    def fit_spe(self,spe_folder,result_name='comptt_gabslog',model='comptt_gabslog',
                run_xspec=0,xspec_comms=[''],**kwargs):


        os.chdir(self.products_path+'/'+spe_folder)
        os.system(f'rm -rf {result_name}.sh')
        specs=glob('*.pi')
        if specs==[]:
            raise Exception('No spectral files found')
        print('Found specta: ',specs)

        xspec_load=nu_xspec.xspec_init(sp1=specs[0],sp2=specs[1],
                                       lowlim=kwargs['lowlim'],uplim=kwargs['uplim'],
                                       folder=result_name)
        run_command(xspec_load,cmd_name=result_name,rewrite=1)


        xspec_model=nu_xspec.xspec_model(model=model,folder=result_name,xspec_comms=xspec_comms)
        run_command(xspec_model,cmd_name=result_name)

        if run_xspec:
            os.system(f'xspec - {result_name}.sh')


    def scan_spe_results(self):
        os.chdir(self.products_path)
        #find all folders with 'spe' in name, scan folders for txt files
        os.system(f'rm -f {self.spe_info_file}')
        for filename in glob('[*spe*,*phase_res*]*/*/*.txt'):
            if len(filename)==0:
                break
            names,values=np.genfromtxt(filename,dtype='str')
            for name,val in zip(names,values):
                self.write_to_obs_info(self.spe_info_file,name,val)


    def _bary_correction_evt(self):
        os.chdir(self.out_path)
        for mode in ['A','B']:
            filename=glob(f'*{mode}01_cl.evt')[0]
            self.barycorr(infile=filename, cmd_path='evt_barycor',barytime='no')

    def _orb_correction_evt(self,mode='A'):
        os.chdir(self.out_path)
        filename=glob(f'*{mode}01_cl.evt_bary')[0]
        q=input(f'Do you want to start doppler correction of {filename}?')
        if q:
            orb_params=doppler_correction.orb_params_v0332
            doppler_correction.correct_times(fitsfile=filename,orb_params=orb_params,time_orig_col='time')

    def orb_correction_lc(self,filename,folder='spe_and_lc',q=1):
        '''
        I use the correction of LC file because the correction of event files
        would take a lot of time
        '''
        os.chdir(self.out_path)
        os.chdir('products')
        os.chdir(folder)

        #filename=glob(f'*{mode}_sr*bary*')[0]
        if q==1:
            q=input(f'Do you want to start doppler correction of {filename}?')
            if q=='y' or q=='yes':
                q=1
            else:
                q=0
        else:
            q=1
        if q:
            orb_params=doppler_correction.orb_params_v0332
            doppler_correction.correct_times(fitsfile=filename,orb_params=orb_params,time_orig_col='time')


    def _make_gti_from_lc(self,folder,mode,
                         period,phase_bins=8,
                         outfolder='gtis/',
                         half_shift=0):

        '''
        This task uses BARYCENTRED AND ORBIRALLY CORRECTED LIGHTCURVES to produce GTI with
        given period and timebins.
        Al long as original lc was barycorrected to different file,
        and this barycentred file was used to produce Orbital correction (see above),
        we juxtapose times from corected LC (TBD time) and original
        LC (times starts with zero). Timezero for the lightcurve is a time of
        the first event in an event file (hence we use time of original lc+ timezero)

        '''
        period=float(period)
        os.chdir(self.out_path)
        os.chdir('products')
        os.chdir(folder)
        gtis_created=glob(f'{outfolder}*{mode}.fits')
        if gtis_created:
            print(gtis_created)
            if input('GTIs for this mode have already been created. Delete folder to rebuild? [y/n]')=='y':
                os.system(f'rm -r {outfolder}')
            else:
                raise Exception('Abort GTI creation')
        create_dir(outfolder)
        filename_orig=glob(f'*{mode}_sr.lc')[0]
        filename_corr=glob(f'*{mode}_sr.bary_lc_orb_corr')[0]

        ff_orig=fits.open(filename_orig)
        time_orig=ff_orig[1].data['time']
        time_orig_zero=ff_orig[1].header['timezero']
        time_orig +=time_orig_zero #no need to add mjdref, because gti is in seconds since ref (see *gti* files in the root directory)
        time_corr=fits.open(filename_corr)[1].data['time']
        if min(np.diff(time_corr))<0:
            raise Exception('corrected times are not monotonous!')


        phase_duration=period/phase_bins
        t0=time_corr[0]
        phases=np.floor((time_corr-t0+half_shift*phase_duration/2)/(period/phase_bins) % phase_bins)+1

        #mindt=np.min(np.diff(time_orig))
        for i in range(1,phase_bins+1):
            ind=start_stop(phases==i)
            gtis=time_orig[ind]

            tstart=gtis[:,0]#-mindt/5
            tstop=gtis[:,1]#+mindt/5


            hdu=fits.PrimaryHDU()
            c1=fits.Column(name='START',format='1D',unit='s',array=tstart)
            c2=fits.Column(name='STOP',format='1D',unit='s',array=tstop)
            cols=fits.ColDefs([c1,c2])
            datahdu=fits.BinTableHDU.from_columns(cols)
            datahdu.name='GTI'
            datahdu.header['START']=len(tstart)
            datahdu.header['STOP']=len(tstop)
            hdulist=fits.HDUList([hdu,datahdu])
            hdulist.writeto(outfolder+f'gti_{i}{mode}.fits')
            hdulist.close()
        print(f'GTIs have been created for module {mode} with period {period} and {phase_bins} bins')


    def _phase_resolved_spectra(self,gtipath='spe_and_lc/gtis',folder='phase_resolved',mode='A'):
        os.chdir(self.out_path)
        os.chdir('products')
        create_dir(folder)
        #outdir=os.path.abspath(folder)
        gtipath=os.path.abspath(gtipath)
        gtis=glob(f'{gtipath}/*{mode}.fits')
        if input(f'GTIs from {gtipath} for module {mode}: \n {[os.path.basename(x) for x in gtis]} \n Is this correct? [y/n]?:  ')=='y':
            os.chdir(folder)
        else:
            raise Exception('Stop phase reloved spectroscopy')

        for ph_num,gtifile in enumerate(gtis,1):
            #self.make_spe(outdir=folder,usrgtifile=gtifile,
            #              stemout=f'phase_resolved_{mode}_bin{ph_num}',
            #              mode=mode,usrgtibarycorr='no')
            self.nuproducts(outdir=folder,usrgtifile=gtifile,
                          stemout=f'phase_resolved_{mode}_bin{ph_num}',
                          mode=mode,usrgtibarycorr='no',
                          barycor=1,grppha=1)
            #self.make_lc(outdir=folder,usrgtifile=gtifile,
            #              stemout=f'phase_resolved_{mode}_bin{ph_num}',
            #              mode=mode,usrgtibarycorr='no')

    def fasebin(self):
        for mode in ['A','B']:
            ser=self.pandas_series()
            per_col=[col for col in ser.index.values if 'period' in col and 'lc'+mode in col and 'err' not in col]
            period_mode=ser[per_col].values[0]
            period_mode='%.4f'% period_mode
            print('PERIOD '+mode, period_mode)
        period_mode=input('Enter period value for phase_resolved spectroscopy: ')
        for mode in ['A','B']:
                self._make_gti_from_lc('spe_and_lc',mode=mode,period=period_mode)
                self._phase_resolved_spectra(mode=mode)



    def fit_ph_res_spe(self,spe_folder='phase_resolved',
                       result_name='comptt_gabslog',model='comptt_gabslog',
                xspec_comms=[''],**kwargs):


        os.chdir(self.products_path+'/'+spe_folder)
        os.system(f'rm -rf {result_name}.sh')
        specs=glob('*.pi')
        nph=int(len(specs)/2)


        if specs==[]:
            raise Exception('No spectral files found')
        print('Found specta: ',specs)


        for phasebin in range(1,nph+1):
            phasebin_specs=glob(f'*bin{phasebin}*.pi')
            xspec_load=nu_xspec.xspec_init(sp1=phasebin_specs[0],sp2=phasebin_specs[1],
                                           lowlim=kwargs['lowlim'],uplim=kwargs['uplim'],
                                           folder=result_name+f'_bin{phasebin}')
            run_command(xspec_load,cmd_name=result_name,rewrite=0)


            xspec_model=nu_xspec.xspec_model(model,result_name+f'_bin{phasebin}',
                                             xspec_comms)
            run_command(xspec_model,cmd_name=result_name)





    def ph_res_param(self,
                       spe_folder='phase_resolved',
                       model='comptt_gabslog_sigma_free',
                       param='DEFAULT',
                       funct=lambda x:x,
                       nph=8,
                       plot=1, ax=None,**plt_kwargs):
        self.scan_spe_results()
        ser=self.pandas_series()


        collist=ser.index.values
        if param=='DEFAULT':

            params=[col for col in collist if 'bin1' in col and model in col and '_lo' not in col and '_hi' not in col]

            params=[par.replace(model+'_bin1_','') for par in params]
            print('Parameters for this model: ',params)
            return params



        def get_parr_array(ser,parname,funct):

            cols_mean=[col for col in ser.index if 'bin' in col and parname in col and '_lo' not in col and '_hi' not in col and model in col]
            mean=ser[cols_mean].values
            mean=funct(mean)


            cols_lo=[col for col in ser.index if 'bin' in col and parname in col and '_lo'  in col and model in col]

            if not cols_lo:
                return mean, np.vstack((mean-mean,mean-mean))

            lo=ser[cols_lo].values
            lo=funct(lo)

            cols_hi=[col for col in ser.index if 'bin' in col and parname in col and '_hi'  in col and model in col]
            hi=ser[cols_hi].values
            hi=funct(hi)


            hi=(hi-mean)
            lo=(mean-lo)


            err=np.vstack((lo,hi))
            if len(mean)!=nph:
                raise Exception(f'Parameter {parname} has less than {nph} values!')

            return mean,err

        phase=np.arange(0,nph)/nph
        try:
            mean,err=get_parr_array(ser,param,funct)
        except:
            raise Exception(f'Error with get par {param}')
        phase=np.concatenate((phase,phase+1))
        mean=np.concatenate((mean,mean))
        err=np.hstack((err,err))

        if plot:
            if ax==None:
                fig,ax=plt.subplots(figsize=(12,4))
                ax.set_title(f'ObsID: {self.ObsID}, model: {model}, parameter: {param}')
            ax.errorbar(phase,mean,yerr=err,drawstyle='steps-mid',**plt_kwargs)

        return phase,mean,err

    def ph_res_results(self,model='comptt_gabslog_sigma_free',rewrite=1):
        '''
        comptt:
        ph_res_params=['exposure', 'MJD', 'chi2_red', 'dof', 'factor1',
'LineE2', 'Sigma3', 'norm4',
'D5', 'Ecycle6', 'sigma7',
'D8', 'Ecycle9', 'sigma10',
'T012', 'kT13', 'taup14', 'norm16', 'factor17',
'eqw_gaussian', 'flux_gabslog_12_79', 'flux_gabslog_4_12',
 'flux_gabslog_4_79', flux_gabslog_7_12, 'flux_gaussian_4_12']

        cutoffpl:
            ['exposure', 'MJD', 'chi2_red', 'dof', 'factor1',
            'LineE2', 'Sigma3', 'norm4',
            'PhoIndex5', 'HighECut6', 'norm7', 'factor8',
            'eqw_gaussian', 'flux_cutoffpl_4_12', flux_cutoffpl_7_12, 'flux_gaussian_4_12']
        '''
        savepath='/Users/s.bykov/work/xray_pulsars/nustar/plots_results/'
        os.chdir(self.out_path+'/products/phase_resolved/')

        matplotlib.rcParams['figure.figsize'] = 6.6, 6.6
        matplotlib.rcParams['figure.subplot.left']=0.1
        matplotlib.rcParams['figure.subplot.bottom']=0.1
        matplotlib.rcParams['figure.subplot.right']=0.9
        matplotlib.rcParams['figure.subplot.top']=0.85

        self.scan_spe_results()
#        ser=self.pandas_series()

        if rewrite:
            fname=self.fasebin_info_file
            os.system(f'> {fname}') #remove content from fasebin_info

        if 'comptt_gabslog' in model:

            fig = plt.figure()
            rows=8
            cols=3
            ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
            ax_flux=ax_eqw.twinx()

            ax_eqw.set_title(self.ObsID+f'\n  model: {model}')
            self.ph_res_param(model=model,param='eqw_gauss',funct=lambda x: 1000*x,
                              ax=ax_eqw,color='r',alpha=1)
            ax_eqw.set_ylabel('Iron line eqw',color='r')

            self.ph_res_param(model=model,param='flux_gabslog_7_12',funct=lambda x: 10**x/1e-8,
                              ax=ax_flux,color='k',alpha=0.6)
            ax_flux.set_ylabel('flux 7-12 ',color='k')

            ax_sigma=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)

            ax_norm=ax_sigma.twinx()

            self.ph_res_param(model=model,param='Sigma3',funct=lambda x: x,
                              ax=ax_sigma,color='b',alpha=0.8)
            ax_sigma.set_ylabel('Iron line Sigma',color='b')

            self.ph_res_param(model=model,param='norm4',funct=lambda x: 1000*x,
                              ax=ax_norm,color='g',alpha=0.6)
            ax_norm.set_ylabel('Iron line norm *1000',color='g')


            ax_chi=plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)
            self.ph_res_param(model=model,param='chi2_red',
                              ax=ax_chi,color='m',alpha=0.9)
            ax_chi.set_ylabel('chi^2/dof',color='m')


            ax_T0 = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=3)
            ax_kT=ax_T0.twinx()

            self.ph_res_param(model=model,param='T012',funct=lambda x: x,
                              ax=ax_T0,color='c',alpha=1)
            ax_T0.set_ylabel('CompTT T0',color='c')

            self.ph_res_param(model=model,param='kT13',funct=lambda x: x,
                              ax=ax_kT,color='m',alpha=1)
            ax_kT.set_ylabel('CompTT kT',color='m')

        if 'cutoffpl' in model:
            fig = plt.figure()
            rows=8
            cols=3
            ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
            ax_flux=ax_eqw.twinx()

            ax_eqw.set_title(self.ObsID+f'\n  model: {model}')
            self.ph_res_param(model=model,param='eqw_gauss',funct=lambda x: 1000*x,
                              ax=ax_eqw,color='r',alpha=1)
            ax_eqw.set_ylabel('Iron line eqw',color='r')

            self.ph_res_param(model=model,param='flux_cutoffpl_7_12',funct=lambda x: 10**x/1e-8,
                              ax=ax_flux,color='k',alpha=0.6)
            ax_flux.set_ylabel('flux 7-12 ',color='k')

            ax_sigma=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)

            ax_norm=ax_sigma.twinx()

            self.ph_res_param(model=model,param='Sigma3',funct=lambda x: x,
                              ax=ax_sigma,color='b',alpha=0.8)
            ax_sigma.set_ylabel('Iron line Sigma',color='b')

            self.ph_res_param(model=model,param='norm4',funct=lambda x: 1000*x,
                              ax=ax_norm,color='g',alpha=0.6)
            ax_norm.set_ylabel('Iron line norm *1000',color='g')


            ax_chi=plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)
            self.ph_res_param(model=model,param='chi2_red',
                              ax=ax_chi,color='m',alpha=0.9)
            ax_chi.set_ylabel('chi^2/dof',color='m')

            #'PhoIndex5', 'HighECut6'
            ax_PI = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=3)
            ax_ecut=ax_PI.twinx()

            self.ph_res_param(model=model,param='PhoIndex5',funct=lambda x: x,
                              ax=ax_PI,color='c',alpha=1)
            ax_PI.set_ylabel('phot. index',color='c')

            self.ph_res_param(model=model,param='HighECut6',funct=lambda x: x,
                              ax=ax_ecut,color='m',alpha=1)
            ax_ecut.set_ylabel('E_cut',color='m')



        # ax = plt.gca()
        # ax.grid(axis='x')
        fig.tight_layout()
        sns.despine(fig,top=1,right=0)
        plt.savefig(savepath+f'ph_res_{self.ObsID}_{model}.png',dpi=500)

        # ax_chi=plt.subplot2grid((rows,cols), (1, 0), rowspan=1, colspan=3)
        # ax_ccf=plt.subplot2grid((rows,cols), (3, 2), rowspan=2, colspan=3)
        # ax_ecutpars=plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3)
        # ax_flux=plt.subplot2grid((rows,cols), (0, 3), rowspan=2, colspan=2)
        # ax_per=ax_flux.twinx()
        # #ax_chi_hist=plt.subplot2grid((rows,cols), (2, 3), rowspan=1, colspan=1)
        # ax_per_find=plt.subplot2grid((rows,cols), (2, 4), rowspan=1, colspan=1)
        # ax_po_and_norm=plt.subplot2grid((rows,cols), (5, 0), rowspan=1, colspan=3)







