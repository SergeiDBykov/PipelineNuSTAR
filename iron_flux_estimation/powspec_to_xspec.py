#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:42:22 2020

@author: s.bykov
"""

from PipelineNuSTAR.core import *

os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/lc_sa')

filename='xspec-powspec/powspec_2031_0.1sec_leahy.qdp'

powspec_data=np.genfromtxt(filename,skip_header=3,usecols=(0,1,2,3))

f,df,P,dP=powspec_data[:,0],powspec_data[:,1],powspec_data[:,2],powspec_data[:,3]
'''
f, df, P, dP

 This should be converted to a text file with the four columns
 being f-df, f+df, 2*P*df, 2*dP*df.
 This file can be used as the input to flx2xsp to create a "spectrum"

'''
xspec_data=np.copy(powspec_data)
xspec_data[:,0]=f-df
xspec_data[:,1]=f+df
xspec_data[:,2]=2*P*df
xspec_data[:,3]=2*dP*df

np.savetxt('xspec-powspec/powspec_2031_xspec.txt',xspec_data,delimiter=' ')