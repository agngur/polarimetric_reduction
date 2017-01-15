from alipy import pysex
import os
from stacks import create_coo_boxes, select_exp
import glob
import numpy as np
from astropy.table import Table, Column, Row
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.utils import random_cmap


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
            try:
                conf_args[name] = value
            except KeyError:
                print('Wrong key for flux measure: {}'.format(name))

    cat = pysex.run(im, im_ref, keepcat=False, sex_command='sextractor',
                    rerun=True,
                    params=['X_IMAGE', 'Y_IMAGE', 'FLUX_BEST', 'FLUXERR_BEST'],
                       conf_args=conf_args)
    
    return cat


def make_plots(im_ref, table):
    


def save_PD_catalog(coo_cats, im_ref, visu=True):
    
    Q = (coo_cats[1]['FLUX_BEST']-coo_cats[0]['FLUX_BEST']) / \
            (coo_cats[1]['FLUX_BEST']+coo_cats[0]['FLUX_BEST'])
    U = (coo_cats[3]['FLUX_BEST']-coo_cats[2]['FLUX_BEST']) / \
            (coo_cats[3]['FLUX_BEST']+coo_cats[2]['FLUX_BEST'])
    PD = np.sqrt(Q**2 + U**2) * 100
    
    table = Table(data=[coo_cats[0]['X_IMAGE'],
                        coo_cats[0]['Y_IMAGE'],
                        Q, U, PD],
                  names=['X_IMAGE', 'Y_IMAGE', 'Q', 'U', 'PD'])
    
    table.sort('PD')
    table_name = im_ref.replace('P1_', '')
    table_name = table_name.replace('.fit', '.cat')
    ascii.write(table, table_name, overwrite=True)
    
    if visu:
        make_plots(im_ref, table)


def make_photometry(work_dir, files_ext='*.fit', ext='.fit',
                    filters=['P1', 'P2', 'P3','P4'],
                    exps=[5, 60],
                    output_dir=None):

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
                    cat = flux_measure(im, im_ref)
                    coo_cats.append(cat)

                save_PD_catalog(coo_cats, im_ref, visu=True) 
