#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 19:52:25 2020

@author: s.bykov
"""


from PipelineNuSTAR.core import *

os.chdir('/Users/s.bykov/work/xray_pulsars/nustar/results/out80102002002/products/lc712_1')


ff1=fits.open('lc712AB_sr.lc_bary_orb_corr')
time=ff1[1].data['time']
rate1=ff1[1].data['rate']

ff2=fits.open('fe_line_0.75.lc_bary_orb_corr')
rate2=ff2[1].data['rate']


CCF_obj_crosscorr=cross_correlation.CrossCorrelation(time, rate1,rate2,circular=0)
CCF_obj_crosscorr.calc_ccf()

fig,ax=plt.subplots(1,sharex='all')
ax.plot(CCF_obj_crosscorr.lag,CCF_obj_crosscorr.ccf,'b:.')
ax.grid()
plt.xlim(-30,30)
plt.show()
