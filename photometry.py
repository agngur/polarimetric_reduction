from alipy import pysex
import os
from stacks import create_coo_boxes, select_exp
import glob
import numpy as np
from astropy.table import Table, Column, Row
from astropy.io import ascii, fits
from astropy.wcs import WCS, utils
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5, ICRS
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch, HistEqStretch, MinMaxInterval, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.utils import random_cmap
from photutils import CircularAperture


def flux_measure(im, im_ref, **kwargs):
    
    conf_args={'VERBOSE_TYPE': 'QUIET',
               'BACKPHOTO_TYPE': 'GLOBAL',
               'DETECT_THRESH': 100.0,
               'ANALYSIS_THRESH': 100.0,
               'GAIN': 1.0,
               'CLEAN': 'Y',
               'SATUR_LEVEL': '60000.0',
               'ASSOC_DATA': '1,2,3',
               'ASSOC_TYPE': 'NEAREST',
               'ASSOCSELEC_TYPE': 'MATCHED'}
    
    if kwargs:
        for name, value in kwargs.items():
            if name in conf_args:
                conf_args[name] = value
            
            """
            try:
                conf_args[name] = value
            except KeyError:
                print('Wrong key for flux measure: {}'.format(name))
            """
            
    cat = pysex.run(im, im_ref, keepcat=False,
                    rerun=True,
                    params=['X_IMAGE', 'Y_IMAGE', 'FLUX_BEST', 'FLUXERR_BEST', 'BACKGROUND', 'THRESHOLD'],
                       conf_args=conf_args)
    
    return cat


def make_plots(im_ref, table):

    plot_name = im_ref.replace('P1_', '')
    plot_name = plot_name.replace('.fit', '.png')
    
    norm = ImageNormalize(stretch=SqrtStretch(), interval=ZScaleInterval())
    fig, ax1 = plt.subplots(1, 1)
    data = fits.getdata(im_ref)
    data_mean, data_std = np.mean(data), np.std(data)
    
    #ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    ax1.imshow(data, interpolation='nearest', cmap='gray', vmin=data_mean-data_std,
               vmax=data_mean+data_std, origin='lower')
    
    
    for star in table:
        ap = CircularAperture([star['X_IMAGE']-1, star['Y_IMAGE']-1], r=4.)
        if star['PD'] < 5:
            color = 'blue'
        elif star['PD'] < 10:
            color = 'yellow'
        else:
            color = 'red'
        
        ap.plot(color=color, lw=0.5, alpha=0.5, ax=ax1)
    
    fig.savefig(plot_name, dpi=300)
    plt.close(fig)


def save_PD_catalog(coo_cats, im_ref, visu=True):
    
    Q = (coo_cats[1]['FLUX_BEST']-coo_cats[0]['FLUX_BEST']) / \
            (coo_cats[1]['FLUX_BEST']+coo_cats[0]['FLUX_BEST'])
    U = (coo_cats[3]['FLUX_BEST']-coo_cats[2]['FLUX_BEST']) / \
            (coo_cats[3]['FLUX_BEST']+coo_cats[2]['FLUX_BEST'])
    PD = np.sqrt(Q**2 + U**2) * 100
    
    #add Q_err i U_err i PD_err
    Q_ERR = np.sqrt((4*((coo_cats[1]['FLUX_BEST'])**2*(coo_cats[0]['FLUXERR_BEST'])**2+ \
                         (coo_cats[0]['FLUX_BEST'])**2*(coo_cats[1]['FLUXERR_BEST'])**2))/ \
                      ((coo_cats[0]['FLUX_BEST']+coo_cats[1]['FLUX_BEST'])**4))
    U_ERR = np.sqrt((4*((coo_cats[3]['FLUX_BEST'])**2*(coo_cats[2]['FLUXERR_BEST'])**2+ \
                         (coo_cats[2]['FLUX_BEST'])**2*(coo_cats[3]['FLUXERR_BEST'])**2))/ \
                      ((coo_cats[2]['FLUX_BEST']+coo_cats[3]['FLUX_BEST'])**4))
    PD_ERR = np.sqrt((Q**2*Q_ERR**2+U**2*U_ERR**2)/(Q**2+U**2))
    
    FLAGS = []
    for i in PD: 
        if i>=100:
            FLAGS.append(int(1))
        else:
            FLAGS.append(int(0))
    
    #statistics to detect/analysis_threshold obtaining
    #max_bgrd = max(coo_cats[0]['FLUXERR_BEST']/coo_cats[0]['FLUX_BEST']*100)
    #print(max_bgrd, '%')
    #max_bgrd = max(coo_cats[1]['FLUXERR_BEST']/coo_cats[1]['FLUX_BEST']*100)
    #print(max_bgrd, '%')
    #max_bgrd = max(coo_cats[2]['FLUXERR_BEST']/coo_cats[2]['FLUX_BEST']*100)
    #print(max_bgrd, '%')
    #max_bgrd = max(coo_cats[3]['FLUXERR_BEST']/coo_cats[3]['FLUX_BEST']*100)
    #print(max_bgrd, '%')
    
    #We can (but esspecially don't want) compute Polarization angle:
    #PA = 1/2 * np.arctan(U/Q)
    #PA_ERR = 1/2*np.sqrt((U**2*Q_ERR**2+Q**2*U_ERR**2)/((Q**2+U**2)**2))
    
    table = Table(data=[coo_cats[0]['X_IMAGE'],
                        coo_cats[0]['Y_IMAGE'], Q, Q_ERR, U, U_ERR, PD, PD_ERR, FLAGS],
                  names=['X_IMAGE', 'Y_IMAGE', 'Q', 'Q_ERR', 'U', 'U_ERR', 'PD', 'PD_ERR', 'FLAGS'])
    
    #FLUX_BEST is the best of FLUX_AUTO and FLUX_ISOCOR
    table_coo = Table(data=[coo_cats[0]['FLUX_BEST'],
                            coo_cats[0]['X_IMAGE'],
                            coo_cats[0]['Y_IMAGE']], 
                      names=['FLUX_BEST','#X_IMAGE', 'Y_IMAGE'])
    
    table.sort('PD')
    table_name = im_ref.replace('P1_', '')
    table_name = table_name.replace('.fit', '.cat')
    ascii.write(table, table_name, overwrite=True)
    
    #file needed to astrometry API
    #table_coo.sort('FLUX_BEST')
    #table_coo.reverse()
    #table_coo.remove_column('FLUX_BEST')
    #table_name = im_ref.replace('P1_', '')
    #table_name = table_name.replace('.fit', '.astrometry')
    #ascii.write(table_coo, table_name, overwrite=True)
    
    if visu:
        make_plots(im_ref, table)

def prepare_catalog(work_dir):
    
    fit_astrometry_files = glob.glob(work_dir+'*P1*'+'*.new')
    for fit_file in fit_astrometry_files:
    
        f = fits.open(fit_file)
        mywcs = WCS(f[0].header)
        f.close()
    
        table_name = fit_file.replace('P1_', '')
        table_file = table_name.replace('.new', '.cat')
        t = Table.read(table_file, format='ascii')
    
        ra_deg = []
        dec_deg = []
        ra_hms = []
        dec_dms = []
        for i in range(len(t)):
            ra, dec = mywcs.all_pix2world([[t['X_IMAGE'][i], t['Y_IMAGE'][i]]], 0)[0]
            c = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
            ra_deg.append(c.ra.deg)
            dec_deg.append(c.dec.deg)
            ra_hms.append(c.ra.to_string(unit=u.hourangle, sep=':', precision=3, pad=True))
            dec_dms.append(c.dec.to_string(sep=':', precision=3, alwayssign=True, pad=True))

        table = Table(data=[t['X_IMAGE'], t['Y_IMAGE'], ra_deg, dec_deg, ra_hms, dec_dms, t['Q'], t['Q_ERR'], t['U'], t['U_ERR'], t['PD'], t['PD_ERR'],t['FLAGS']], 
                      names=['X_IMAGE', 'Y_IMAGE','RA_deg', 'DEC_deg', 'RA', 'DEC', 'Q', 'Q_ERR', 'U', 'U_ERR', 'PD', 'PD_ERR', 'FLAGS'])    

        table.write(table_file.replace('.cat', '.catalog'), format='ascii', overwrite=True)
        
def prepare_catalog(work_dir):
    
    fit_astrometry_files = glob.glob(work_dir+'*P1*'+'*.new')
    for fit_file in fit_astrometry_files:
    
        f = fits.open(fit_file)
        mywcs = WCS(f[0].header)
        f.close()
    
        table_name = fit_file.replace('P1_', '')
        table_file = table_name.replace('.new', '.cat')
        t = Table.read(table_file, format='ascii')
    
        ra_deg = []
        dec_deg = []
        ra_hms = []
        dec_dms = []
        for i in range(len(t)):
            ra, dec = mywcs.all_pix2world([[t['X_IMAGE'][i], t['Y_IMAGE'][i]]], 0)[0]
            c = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
            ra_deg.append(c.ra.deg)
            dec_deg.append(c.dec.deg)
            ra_hms.append(c.ra.to_string(unit=u.hourangle, sep=':', precision=3, pad=True))
            dec_dms.append(c.dec.to_string(sep=':', precision=3, alwayssign=True, pad=True))

        table = Table(data=[t['X_IMAGE'], t['Y_IMAGE'], ra_deg, dec_deg, ra_hms, dec_dms, t['Q'], t['Q_ERR'], t['U'], t['U_ERR'], t['PD'], t['PD_ERR'],t['FLAGS']], 
                      names=['X_IMAGE', 'Y_IMAGE','RA_deg', 'DEC_deg', 'RA', 'DEC', 'Q', 'Q_ERR', 'U', 'U_ERR', 'PD', 'PD_ERR', 'FLAGS'])    

        table.write(table_file.replace('.cat', '.catalog'), format='ascii', overwrite=True)


def make_photometry(work_dir, files_ext='*.fit', ext='.fit',
                    filters=['P1', 'P2', 'P3','P4'],
                    exps=[5, 60],
                    output_dir=None, 
                    **kwargs):

    if output_dir == None:
        output_dir = os.path.join(work_dir, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    time = []
    images = sorted(glob.glob(os.path.join(work_dir, files_ext)))
    for exp in exps:
        selected_exp = select_exp(images, exp, ext, prefix='_')
        coo_boxes = create_coo_boxes(selected_exp)
        for coo_box in coo_boxes:
            coo_box.sort(key=lambda name: os.path.basename(name).split('_')[-2])
            im_ref = coo_box[0]
            f = fits.open(coo_box[0], mode='update')
            time.append(f[0].header['JD'])
            if len(coo_box) < 4:
                print('coo box too short')
            else:
                coo_cats = []
                for im in coo_box:
                    cat = flux_measure(im, im_ref, **kwargs)
                    coo_cats.append(cat)
                coo_cats.append(time)
                print(len(coo_cats[1]), len(time))
                print(coo_cats)

                save_PD_catalog(coo_cats, im_ref, visu=True) 
 
