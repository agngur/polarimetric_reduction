from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
import alipy
import numpy as np
import os
import glob
import datetime as dt
import processimages as pim
import logging


def start_log(main_dir):
    logging.basicConfig(filename=os.path.join(main_dir, 'prepare_stack.log'),
                        format='%(asctime)s %(message)s', level=logging.INFO)

    
def create_im_list(main_dir, start_date, end_date, ext):
    container = []
    im_amount = 0
    for root, dirs, files in os.walk(main_dir, topdown=False):
        for name in dirs:
            folder_root_name = os.path.join(root, name)
            # format w zaleznosci od systemu
            # folder_name = folder_root_name.split('\\')[-1]
            folder_name = folder_root_name.split('/')[-1]
            try:
                f_date = dt.datetime.strptime(folder_name, '%Y-%m-%d')
                if (f_date >= start_date and f_date <= end_date):
                    folder_content = sorted(glob.glob(
                        os.path.join(folder_root_name, ext)))
                    if len(folder_content) > 0:
                        container.append([folder_root_name, folder_content])
                        im_amount += len(folder_content)
            except ValueError:
                pass
            
    logging.info('Found {:d} reducted images'.format(im_amount))       
    return container


def select_exp(pack, exp, files_ext, prefix):
    express = ''.join([prefix, str(exp), files_ext])

    return list((filter(lambda im: express in im, pack)))


def select_filter(pack, fil):
    
    return list(filter(lambda im: fil in im, pack))


def create_coo_boxes(pack, radius=4, unit='deg' , method='name'):

    images_fields = []
    grouped_fields = []

    if method == 'name':
        for im in pack:
            im_name = os.path.basename(im)
            im_field = im_name.split('_')[1]
            images_fields.append([im, im_field])

        for im_field in images_fields:
            group = [f[0] for f in images_fields if f[1] == im_field[1]]
            if len(group) > 0: grouped_fields.append(group)
            images_fields = [f for f in images_fields if f[0] not in group]
    
        return grouped_fields
    
    elif method == 'coo':
    
        radius = radius * getattr(u, unit)
    
        for im in pack:
            im_name = os.path.basename(im)
            # sensitive to name change
            im_field = im_name.split('_')[1]
        
            if '+' in im_field:
                im_ra, im_dec = im_field.split('+')
            elif '-' in im_field:
                im_ra, im_dec = im_field.split('-')
            else:
                logging.error('Something wrong with coordinates format in im name')
                break
            
            # one line, slightly complicated solution for change HHMMSS to HH:MM:SS 
            im_ra = ':'.join(list(map(''.join, zip(*[iter(im_ra)]*2))))
            im_dec = ':'.join(list(map(''.join, zip(*[iter(im_dec)]*2))))
            im_field = SkyCoord(ra=im_ra, dec=im_dec, unit=(u.hourangle, u.deg))
            images_fields.append([im, im_field])
        
        for im_field in images_fields:
            group = [f[0] for f in images_fields if f[1].separation(im_field[1]) < radius]
            if len(group) > 0: grouped_fields.append(group)
            images_fields = [f for f in images_fields if f[0] not in group]
        
        return grouped_fields
    
    else:
        logging.error('Wrong method name')
        brake
        

def clean_unaligned(raw_aligned_filter_packs):

    suffixes = []
    for i, pack in enumerate(raw_aligned_filter_packs):
        suffixes.append('_'.join(pack[0].split('_')[3:]))
        raw_aligned_filter_packs[i] = set('_'.join(x.split('_')[:3]) for x in pack)

    intersection = set.intersection(*raw_aligned_filter_packs)

    aligned_filter_packs = []
    for suffix in suffixes:
        aligned_filter_packs.append(['_'.join((x, suffix)) for x in intersection])

    return aligned_filter_packs


def prepare_stack(main_dir, save_dir, hdr_keys, start_date, end_date, files_ext='_red.fit',
                  ext='*red.fit', exps=[5], 
                  filters=['P1', 'P2', 'P3', 'P4'], 
                  radius=4, unit='deg', logs=True):
    print('hello')
    try:
        os.mkdir(save_dir)

    except FileExistsError:
        pass
    
    if logs:
        start_log(main_dir)
        logging.info('PREPARE STACK START: '
                     'main dir={}'
                     'save dir={}'.format(main_dir, save_dir))
        
    container = create_im_list(main_dir, start_date, 
                               end_date, ext=ext)
    

    for pack_date, pack in container:
        for exp in exps:
            selected_exp = select_exp(pack, exp, files_ext, prefix='-') 
            logging.info('Selected {:d} images with exp {}'.format(len(selected_exp), exp))
            if selected_exp:
                coo_boxes = create_coo_boxes(selected_exp, radius, unit, method='name')
                logging.info('Selected {:d} coo groups'.format(len(coo_boxes)))
                for coo_box in coo_boxes:
                    raw_aligned_filter_packs = []
                    stacked_images = []
                    for filter_name in filters:
                        filter_pack = select_filter(coo_box, filter_name)
                        logging.info('Selected {:d} images with filter {}'.format(len(filter_pack), 
                                                                                  filter_name))
                        filter_pack = pim.align_images(filter_pack)
                        raw_aligned_filter_packs.append(filter_pack)
                        logging.info('Pack align done')
                        
                    aligned_filter_packs = clean_unaligned(raw_aligned_filter_packs)
                    for filter_name, filter_pack in zip(filters, aligned_filter_packs):
                        stacked_image = pim.make_stack(filter_pack, save_dir,
                                                           exp, filter_name, hdr_keys)
                        logging.info('Pack stack done')
                        stacked_images.append(stacked_image)

                    pim.align_images(stacked_images)
                    logging.info('Masters align done')

    
