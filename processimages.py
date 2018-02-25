# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Fri Jan  6 21:03:39 2017
@author: Agnieszka & Michal
"""

import gzip
import glob
import astropy.io.fits as fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
import subprocess as sub
from ccdproc import CCDData, Combiner, subtract_bias, subtract_dark, flat_correct
import numpy as np
import os
from shutil import copyfile
import alipy
import ccdproc
import time 

from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [20, 20]

#-----------------------------------------------------#
# solve_filed param
time_limit = 15     # time limit before giving up
scale_low = 32      # arcsec
scale_high = 34     # arcsec
solve_radius = 5    # deg
#solve_depth = '40,80,100,160,250'

# clear
files_to_rm = ['*.axy', '*.corr', '*.xyls', '*.match', '*.rdls', '*.solved', '*.wcs', '*objs.png']

#-----------------------------------------------------#


def uint_all(im_dir):
    all_fit_files = glob.glob(os.path.join(im_dir, '*.fit'))
    for file in all_fit_files:
        f = fits.open(file, mode='update')
        f[0].data = f[0].data.astype(np.uint16)
        f.flush()
        f.close()
        

#im_dir = 'C:(...)'
def create_MasBias(im_dir):
 
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*bias.fit')))
    #print(images)

    for im in images:
        CCD_data_table.append(CCDData.read(im, unit='adu'))
    
    combiner = Combiner(CCD_data_table, dtype='int16')
    median = combiner.median_combine()

    CCDData.write(median, os.path.join(im_dir, 'MBias.fit'),  hdu_mask=None, hdu_uncertainty=None, clobber=True)
     

def create_MasDark(im_dir, dark_exp, m_bias_name='MBias.fit'):
    
    m_bias =  CCDData.read(m_bias_name)
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*dark'+str(dark_exp)+'.fit')))

    for im in images:
        dark_data = CCDData.read(im, unit='adu')
        dark_data = subtract_bias(dark_data, m_bias)
        CCD_data_table.append(dark_data)
        
    combiner = Combiner(CCD_data_table, dtype='int16')
    #median = combiner.median_combine()
    median = combiner.sigma_clipping()
    
    CCDData.write(median, os.path.join(im_dir, 'MDark'+str(dark_exp)+'.fit'), hdu_mask=None, 
                  hdu_uncertainty=None, clobber=True)


def create_MasFlat(im_dir, dark_exp, m_dark_name, m_bias_name='MBias.fit'):
    
    m_bias =  CCDData.read(m_bias_name, unit='adu')
    m_dark =  CCDData.read(m_dark_name, unit='adu')
    CCD_data_table = []
    images = sorted(glob.glob(os.path.join(im_dir, '*flatR.fit')))
    median_table = []
    for im in images:
        flat_data = CCDData.read(im, unit='adu')
        median_table.append(np.median(flat_data.data))
        flat_data = subtract_bias(flat_data, m_bias)
        flat_data = subtract_dark(flat_data, m_dark, 
                                  dark_exposure=dark_exp*u.second, 
                                  data_exposure=float(
                                        flat_data.header['EXPTIME'])*u.second)
        flat_data.data = flat_data.data / np.median(flat_data.data) 
        CCD_data_table.append(flat_data)
         
    combiner = Combiner(CCD_data_table, dtype=np.float64)
    median = combiner.median_combine()
    scale = np.median(median_table)
    print(scale)
    print(np.median(median.data))
    median.data = median.data*scale
    median.data = median.data.astype(np.uint16)
    
    CCDData.write(median, os.path.join(im_dir, 'MFlat.fit'), hdu_mask=None, hdu_uncertainty=None, clobber=True)
    

def reduction_data(im_dir, exp_time, dark_exp, m_dark_name, 
                   m_flat_name='MFlat.fit', m_bias_name='MBias.fit'):

    m_bias =  CCDData.read(m_bias_name, unit='adu')
    m_dark =  CCDData.read(m_dark_name, unit='adu')
    m_flat =  CCDData.read(m_flat_name, unit='adu')
    m_flat.data = m_flat.data / np.max(m_flat.data)
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
        CCDData.write(ccd_data, im.replace('.fit', '_red.fit'), hdu_mask=None, 
                      hdu_uncertainty=None, clobber=True)
        '''
        f = fits.open(im.replace('.fit', '_red.fit'), mode='update')
        f[0].header = fits.getheader(im)
        f.flush()
        '''

def make_stack(images, save_dir, exp, filter_name, hdr_keys, ref_im_num=1):
    
    print('images', images)
    # CCD_data_table = list(map(CCDData.read(unit='adu'), images))
    CCD_data_table = [CCDData.read(im, unit='adu') for im in images]
    combiner = Combiner(CCD_data_table, dtype='float32')
    # combiner.minmax_clipping(max_clip=50000)
    median = combiner.median_combine()
    
    ref_image_name = os.path.basename(images[ref_im_num]).split('_')
    ref_image_name = '_'.join(ref_image_name[0:2])
    stack_name = '_'.join([ref_image_name, filter_name, str(exp)]) + '.fit'
    master_hdr = create_stack_hdr(images, hdr_keys)
    print(stack_name)
    CCDData.write(median, os.path.join(save_dir, stack_name),
                  hdu_mask=None, hdu_uncertainty=None, clobber=True)
    f = fits.open(os.path.join(save_dir, stack_name), mode='update')
    f[0].header = master_hdr
    f.flush()
    
    return os.path.join(save_dir, stack_name)
                  
    
def create_stack_hdr(filter_pack, hdr_keys):
    
    hdr_keys = np.genfromtxt(hdr_keys, dtype=str)
    h1 = fits.getheader(filter_pack[0])
    h2 = fits.getheader(filter_pack[-1])
    master = fits.Header()

    for key in hdr_keys:
        if key == 'DATE-OBS':
            master['OBSSTART'] = h1[key]
            master['OBSEND'] = h2[key]
            h1_time = Time(h1[key])
            h2_time = Time(h2[key])
            mid_time = h1_time + (h2_time - h1_time) / 2. 
            master[key] = mid_time.value
            
        elif key == ('OBJCTRA' or 'OBJCTDEC'):
            im_coo = SkyCoord(ra=h1['OBJCTRA'], dec=h1['OBJCTDEC'], unit=(u.hourangle, u.deg))
            master['OBJCTRA'] = im_coo.ra.deg
            master['OBJCTDEC'] = im_coo.dec.deg

        else:      
            try:
                float(h1[key])
                master[key[:3]+'_sta'] = h1[key]
                master[key[:3]+'_end'] = h2[key]
                print(h2[key])
                master[key] = (h2[key] + h1[key]) / 2.
            except ValueError:
                master[key] = h1[key]
    master['IMNUM'] = len(filter_pack)

    return master
    
    
def align_images(images, ref_im_num=0, save_dir=None, overwrite=True):
    
    ref_image = images[ref_im_num]
    tempfile = os.path.join(os.path.dirname(ref_image), ".temp")
    copyfile(ref_image, tempfile)
    ref_image = tempfile
   
    identifications = alipy.ident.run(ref_image, images, visu=False, sexkeepcat=True, verbose=False)  
    #sex_command='sextractor'

    for id in identifications: 
        if id.ok == True: 
            print("{} : {}, flux ratio {}".format(id.ukn.name, id.trans, id.medfluxratio))
        else:
            images.remove(id.ukn.filepath)
            print("{} : no transformation found !".format(id.ukn.name))

    outputshape = alipy.align.shape(ref_image)

    for id in identifications:
        if id.ok == True:
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=False,
                                    overwrite=True)

    return images            
            
def do_astrometry(cat_dir, limit=100):
    cat_files = glob.glob(os.path.join(cat_dir,'*.pysexcat'))
    for cat_file in cat_files:
        cat = np.genfromtxt(cat_file, usecols=[0, 1, 2])
        cat = cat[cat[:,2].argsort()[::-1]]
        np.savetxt(cat_file, cat[:limit,:2], fmt="%.4f")

    fit_files = glob.glob(os.path.join(cat_dir,'*.fit'))
    for fit_file in fit_files:
        fits_coo = open_file(fit_file)
        solve_field(fits_coo, fit_file)
        
    clear(cat_dir)

def open_file(fit_file):
    hdr = fits.getheader(fit_file)
    #print('fits coo:', str(hdr['OBJCTRA'])+" "+str(hdr['OBJCTDEC']))
    fits_coo = SkyCoord(str(hdr['OBJCTRA'])+" "+str(hdr['OBJCTDEC']),
                        'icrs', unit=(u.hour, u.deg))
    return fits_coo

def solve_field(fits_coo, fit_file):
    solve_field_command = ['solve-field',
                           '--ra', '%s' % (fits_coo.ra.deg),
                           '--dec', '%s' % (fits_coo.dec.deg),
                           '--radius', '%1.1f' % solve_radius,
                           #'--depth', solve_depth,
                           '--cpulimit', '%f' % time_limit,
                           '--scale-units', 'arcsecperpix',
                           '--scale-low', '%.5f' % scale_low,
                           '--scale-high', '%.5f' % scale_high,
                           '--overwrite',
                           '--no-verify',
                           '--no-plots',
                           str(fit_file)]
    sub.Popen(solve_field_command, stdout=sub.PIPE,
              stderr=sub.PIPE).communicate()
    if os.path.exists(str(fit_file).split(".")[0]+'.new'):
        return True
    else:
        print('solve error')
        return False

def clear(cat_dir):
    print('cleaning.....')
    for i in files_to_rm:
        files = glob.glob(os.path.join(cat_dir, i))
        print(os.path.join(cat_dir, i))
        for j in files:
            os.remove(j)
    print("Files removed!")

    
#create_MasDark(im_dir, 60)
#create_MasDark(im_dir, 5)        
#create_MasFlat(im_dir, 5, 'MDark5.fit')  

#reduction_data(im_dir, 60, 60, 'MDark60.fit')
#reduction_data(im_dir, 5, 5, 'MDark5.fit')

#solve_field(open_file(fits_file), fits_file)
#fit_file = '/home/agngur/Dokumenty/survey/test/2017-06-09/stacked_fieldss/2457914_154312+853514_P1_5.fit'
#solve_field(open_file(fit_file), fit_file)
