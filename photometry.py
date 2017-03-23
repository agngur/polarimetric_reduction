from alipy import pysex
import os
from stacks import create_coo_boxes, select_exp
import glob
import numpy as np
from astropy.table import Table, Column, Row
from astropy.io import ascii, fits
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
                    params=['X_IMAGE', 'Y_IMAGE', 'FLUX_BEST', 'FLUXERR_BEST'],
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
    
    table = Table(data=[coo_cats[0]['X_IMAGE'],
                        coo_cats[0]['Y_IMAGE'], Q, U, PD],
                  names=['X_IMAGE', 'Y_IMAGE', 'Q', 'U', 'PD'])
    
    table_coo = Table(data=[coo_cats[0]['FLUX_BEST'],
                            coo_cats[0]['X_IMAGE'],
                            coo_cats[0]['Y_IMAGE']], 
                      names=['FLUX_BEST','#X_IMAGE', 'Y_IMAGE'])
    
    table.sort('PD')
    table_name = im_ref.replace('P1_', '')
    table_name = table_name.replace('.fit', '.cat')
    ascii.write(table, table_name, overwrite=True)
    
    
    table_coo.sort('FLUX_BEST')
    table_coo.reverse()
    table_coo.remove_column('FLUX_BEST')
    table_name = im_ref.replace('P1_', '')
    table_name = table_name.replace('.fit', '.astrometry')
    
    ascii.write(table_coo, table_name, overwrite=True)
    
    if visu:
        make_plots(im_ref, table)

        
        

def make_photometry(work_dir, files_ext='*.fit', ext='.fit',
                    filters=['P1', 'P2', 'P3','P4'],
                    exps=[5, 60],
                    output_dir=None, 
                    **kwargs):

    if output_dir == None:
        output_dir = os.path.join(work_dir, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    images = sorted(glob.glob(os.path.join(work_dir, files_ext)))
    for exp in exps:
        selected_exp = select_exp(images, exp, ext, prefix='_')
        coo_boxes = create_coo_boxes(selected_exp)
        for coo_box in coo_boxes:
            coo_box.sort(key=lambda name: os.path.basename(name).split('_')[-2])
            im_ref = coo_box[0]
            if len(coo_box) < 4:
                print('coo box too short')
            else:
                coo_cats = []
                for im in coo_box:
                    cat = flux_measure(im, im_ref, **kwargs)
                    coo_cats.append(cat)

                save_PD_catalog(coo_cats, im_ref, visu=True) 
