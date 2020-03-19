#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 13:52:17 2020

@author: s.bykov
"""

#def argslist(args):
#    if not args:
#        return '\n'
#    return ''.join(list(arg+'\n ' for arg in args))

def xspec_init(folder,sp1,sp2,lowlim=4,uplim=79):

    return f'''
#xspec

@/Users/s.bykov/work/my_scripts/xspec/xspec_utils.tcl
data 1:1 {sp1}
data 2:2 {sp2}

query yes
ign *: **-{lowlim}. {uplim}.-**
ign bad
setpl en

rm -rf {folder}/ ; mkdir {folder}/ ;  cd {folder}/

'''

def xspec_model(model,folder,xspec_comms):
    if len(xspec_comms)==1:
        xspec_comms=[xspec_comms[0],'']
    if model=='po':
        load_model=f'''
@/Users/s.bykov/work/xray_pulsars/nustar/python_nu_pipeline/xspec_scripts/po_dim.xcm
'''
        load_model+=''.join(list('\n'+arg+'\n ' for arg in xspec_comms))
        load_model+=f'''
renorm
fit 1000

parallel error 6


error 2,3,4,5,6

writetext {folder}.txt 0 50 {folder}
save model mymodel.xcm


calc_flux 3 4 12 -9  {folder}
calc_flux 2 4 12 -10 {folder}

calc_eqw 2 68 {folder}

cd ../

'''

    if model=='cutoffpl':
        load_model=f'''
@/Users/s.bykov/work/xray_pulsars/nustar/python_nu_pipeline/xspec_scripts/cutoffpl_bright.xcm
'''
        load_model+=''.join(list('\n'+arg+'\n ' for arg in xspec_comms))
        load_model+=f'''
renorm
fit 1000

parallel error 6


error 2,3,4,5,6,7

writetext {folder}.txt 0 50 {folder}

save model mymodel.xcm

calc_eqw 2 68 {folder}


calc_flux 3 4 12 -8  {folder}

calc_flux 3 7 12 -8  {folder}

calc_flux 2 4 12 -10 {folder}

cd ../

'''

    if model=='comptt_gabslog':
        load_model=f'''
@/Users/s.bykov/work/xray_pulsars/nustar/python_nu_pipeline/xspec_scripts/comptt_gabs_bright.xcm
'''
        load_model+=''.join(list('\n'+arg+'\n ' for arg in xspec_comms))
        load_model+=f'''
renorm
fit 1000

parallel error 15


error 2,3,4,5,6,7,12,13,14,16,17 #,8,9,10

writetext {folder}.txt 0 50 {folder}

save model mymodel.xcm


calc_eqw 2 68 {folder}


calc_flux 3 4 12 -8 {folder}

calc_flux 3 7 12 -8 {folder}

calc_flux 3 12 79 -8 {folder}

calc_flux 2 4 12 -10 {folder}


cd ../

'''



    if model=='comptt_gabslog_no_gauss':
        load_model=f'''
@/Users/s.bykov/work/xray_pulsars/nustar/python_nu_pipeline/xspec_scripts/comptt_gabs_no_gauss.xcm
'''
        load_model+=''.join(list('\n'+arg+'\n ' for arg in xspec_comms))
        load_model+=f'''
renorm
fit 1000

parallel error 15


error 2,3,4,5,6,7,9,10,12,13,14

writetext {folder}.txt 0 50 {folder}

save model mymodel.xcm


calc_flux 2 12 79 -8 {folder}

calc_flux 2 4 79 -8 {folder}

calc_flux 2 4 12 -9 {folder}

#cd ../

'''
    return load_model#+'\n exit'



