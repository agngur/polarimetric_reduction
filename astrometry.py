from astrometry_client import Client
from astropy.io import fits
import subprocess as sub
import time
import numpy as np
import os
import glob

def do_astrometry(cat_dir, limit=100):
    cat_files = glob.glob(cat_dir+'/'+'*.pysexcat')
    #for cat_file in cat_files:
    #    cat = np.genfromtxt(cat_file, usecols=[0, 1, 2])
    #    cat = cat[cat[:,2].argsort()[::-1]]
    #    np.savetxt(cat_file, cat[:limit,:2], fmt="%.4f")
        
    fit_files = [f.replace('.pysexcat','.fit') for f in cat_files]
    for fit, cat in zip(fit_files, cat_files):
        make_astrometry(fit, cat)

def make_astrometry(file_name, cat):
    server = 'http://nova.astrometry.net/api/'
    apikey = 'abonxcbacdbfjbce'
    client = Client(apiurl=server)
    client.login(apikey)
    #FIXME
    # im_ref = [im for im in images if 'P1' in im]
    # cat = im_ref[0].replace('_P1', '')
    # cat = cat.replace('.fit', '.astrometry')
    wcs_file = cat.replace('.pysexcat', '.wcs')
    print('Cat ', cat)
    print('WCS: ', wcs_file)
    time.sleep(10)
    cmd = 'python /home/agngur/polarimetric_reduction/polarimetric_reduction/astrometry_client.py\
     --server {server} --apikey {apikey} -u {cat} --wcs={wcs_file}'.format(server=server, apikey=apikey,
                                             cat=cat, wcs_file=wcs_file)
        
    print(cmd)
    #os.system(cmd)
    print('kkkk')
    sub.call(cmd.split(' '))
    #print('sdfdfsdf')
    wcs = fits.getheader(wcs_file)
    for im in [file_name]:
        f = fits.open(im, mode='update')
        f[0].header += wcs
        f.flush()

cat_dir = '/home/agngur/Dokumenty/survey/test/2017-06-09/stacked_fieldss'
do_astrometry(cat_dir)
#make_astrometry('/home/agngur/Dokumenty/survey/test/2017-06-09/stacked_fieldss/2457914_154312+853514_P1_5.fit', 
#                '/home/agngur/Dokumenty/survey/test/2017-06-09/stacked_fieldss/2457914_154312+853514_P1_5.pysexcat')