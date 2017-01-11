# -*- coding: utf-8 -*-
#! usr/bin/python
"""
Created on Fri Jan  6 21:03:39 2017
@author: Agnieszka
"""

import gzip
import glob
from astropy import units as u
from ccdproc import CCDData, Combiner, subtract_bias, subtract_dark, flat_correct
import numpy as np
import os
import ccdproc

#im_dir = 'C:/Users/Agnieszka/Documents/UP/Praca_lic/MARAT/2016-12-21'
def create_MasBias(im_dir):
 
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*bias.fit')))
    #print(images)

    for im in images:
        CCD_data_table.append(CCDData.read(im, unit='adu'))
    
    combiner = Combiner(CCD_data_table, dtype='int16')
    median = combiner.median_combine()

    CCDData.write(median, os.path.join(im_dir, '!MBias.fit'),  hdu_mask=None, hdu_uncertainty=None, clobber=True)
     

def create_MasDark(im_dir, dark_exp, m_bias_name='!MBias.fit'):
    
    m_bias =  CCDData.read(m_bias_name)
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*dark'+str(dark_exp)+'.fit')))

    for im in images:
        dark_data = CCDData.read(im, unit='adu')
        dark_data = subtract_bias(dark_data, m_bias)
        CCD_data_table.append(dark_data)
        
    combiner = Combiner(CCD_data_table, dtype='int16')
    median = combiner.median_combine()
    
    CCDData.write(median, os.path.join(im_dir, '!MDark'+str(dark_exp)+'.fit'), hdu_mask=None, 
                  hdu_uncertainty=None, clobber=True)


def create_MasFlat(im_dir, dark_exp, m_dark_name, m_bias_name='!MBias.fit'):
    
    m_bias =  CCDData.read(m_bias_name, unit='adu')
    m_dark =  CCDData.read(m_dark_name, unit='adu')
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*flatR.fit')))
    #print(images)
    for im in images:
        flat_data = CCDData.read(im, unit='adu')
        flat_data = subtract_bias(flat_data, m_bias)
        flat_data = subtract_dark(flat_data, m_dark, 
                                  dark_exposure=dark_exp*u.second, 
                                  data_exposure=float(
                                        flat_data.header['EXPTIME'])*u.second)
        CCD_data_table.append(flat_data)
        
    combiner = Combiner(CCD_data_table, dtype='int16')
    median = combiner.average_combine()
    
    CCDData.write(median, os.path.join(im_dir, '!MFlat.fit'), hdu_mask=None, hdu_uncertainty=None, clobber=True)


def reduction_data(im_dir, exp_time, dark_exp, m_dark_name, 
                   m_flat_name='!MFlat.fit', m_bias_name='!MBias.fit'):

    m_bias =  CCDData.read(m_bias_name, unit='adu')
    m_dark =  CCDData.read(m_dark_name, unit='adu')
    m_flat =  CCDData.read(m_flat_name, unit='adu')
    #CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, 
                                           '*P*'+'-'+str(exp_time)+'.fit')))

    for im in images:
        ccd_data = CCDData.read(im, unit='adu')
        ccd_data = subtract_bias(ccd_data, m_bias)
        ccd_data = subtract_dark(ccd_data, m_dark, 
                                 dark_exposure=dark_exp*u.second, 
                                 data_exposure=float(
                                        ccd_data.header['EXPTIME'])*u.second)
        ccd_data = flat_correct(ccd_data, m_flat)
        ccd_data.data = ccd_data.data.astype('float32')
        CCDData.write(ccd_data, im.replace('.fit', '_red.fit'), hdu_mask=None, hdu_uncertainty=None, clobber=True)
        

#create_MasDark(im_dir, 60)
#create_MasDark(im_dir, 5)        
#create_MasFlat(im_dir, 5, '!MDark5.fit')  

#reduction_data(im_dir, 60, 60, '!MDark60.fit')
#reduction_data(im_dir, 5, 5, '!MDark5.fit')

