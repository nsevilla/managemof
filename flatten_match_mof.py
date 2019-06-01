#!/usr/bin/env python
'''
flatten_match_mof.py
This script will  create a file for upload containing the flattened data in the files from the mof_dir directory (the MOF files) but only for those objects included in the corresponding id file in coaddid_dir.
Usage: python flatten_match_mof.py -d [mof_dir] -i [coaddid_dir] -o [out_dir]
Author: Nacho Sevilla (nsevilla@gmail.com) with some code from Erin Sheldon
'''
import numpy as np
from astropy.io import fits
import os
import sys
from optparse import OptionParser

bands = ['g','r','i','z']
onedim_cols = ['id','number','fofid','coadd_run','flags','time_last_fit','box_size','obj_flags','mask_frac','psfrec_T','cm_flags','cm_T','cm_T_err','cm_T_s2n','cm_weight','cm_max_flags','cm_max_T','cm_max_T_err','cm_max_T_s2n','cm_s2n_w','cm_chi2per','cm_dof','cm_flags_r','cm_s2n_r','cm_T_r','cm_psf_T_r','cm_fracdev','cm_fracdev_noclip','cm_fracdec_err','cm_TdByTe','cm_TdByTe_noclip','cm_mof_flags','cm_mof_num_itr','fofind','ra','dec']
bidim_cols = ['cm_pars_cov','cm_max_pars_cov','cm_g_cov']
band_cols = ['nimage_tot','psf_flags','psf_flux','psf_flux_err','psf_flux_s2n','psf_mag','nimage_use','cm_flux_cov','cm_flux','cm_flux_s2n','cm_mag','cm_logsb','cm_max_flux_cov','cm_max_flux','cm_max_flux_s2n','cm_max_mag','cm_max_logsb']
bidim_band_cols = ['cm_flux_cov','cm_max_flux_cov','cm_pars_cov','cm_max_pars_cov','cm_g_cov']
multi_cols = ['psfrec_g','cm_pars','cm_g','cm_max_pars', 'cm_mof_abs_diff','cm_mof_frac_diff','cm_mof_err_diff']

def get_coldefs(descr, defs={}, bands=None, band_cols=None, bidim_cols=None):
    """
    Convert a numpy descriptor to a set of oracle 
    column definitions

    array columns are converted to name_{dim1}_{dim2}...{dimn}

    parameters
    ----------
    descr: numpy type descriptor
        E.g. arr.dtype.descr
    defs: dict,optional
        A dict returning a list of field defs. It is keyed by field names from
        the array.  This can be used to over-ride the defaults, e.g. to use a
        different name or to over-ride conversions for arrays.
    """

    if defs is None:
        defs={}

    alldefs=[]
    def_template='%s not null'

    for d in descr:
        name=d[0]
        ot=get_oracle_type(d[1])

        if name in defs:
            alldefs += defs[name]
        elif len(d) == 2:
            # this is a scalar column... easy!
            defi=def_template % ot
            alldefs.append( (name,defi,d[1]) )
        else:
            dims=d[2]
            if not isinstance(dims,tuple):
                dims=(dims,)

            if (bands is not None
                    and band_cols is not None
                    and name in band_cols):
                names=get_band_arr_colnames(name,dims,bands,bidim_cols)
            else:
                names=get_arr_colnames(name,dims)
            
            for n in names:
                defi=def_template % (ot)
                alldefs.append( (n,defi,d[1]) )
    #exit()
    return alldefs

def get_arr_colnames(name, dims):
    """
    Get db names for an array, naming 
        name_{num1}_{num2}...
    """
    ndim=len(dims)
    if ndim==1:
        names=get_arr1_colnames(name,dims)
    elif ndim==2:
        names=get_arr2_colnames(name,dims)
    else:
        raise ValueError("only support 1 and 2 d arrays")

    return names

def get_arr1_colnames(name, dims):
    """
    Get db names for an array, naming 
        name_{num}
    """
    names=[]
    for n in xrange(1,dims[0]+1):
        names.append( '%s_%d' % (name,n) )

    return names

def get_arr2_colnames(name, dims):
    """
    Get db names for an array, naming 
        name_{num1}_{num2}
    """
    names=[]
    for n1 in xrange(1,dims[0]+1):
        for n2 in xrange(1,dims[1]+1):
            names.append( '%s_%d_%d' % (name,n1,n2) )

    return names

def get_band_arr_colnames(name, dims, bands, bidim_cols):
    """
    Get db names for an array, naming 
        name_{num1}_{num2}...
    """
    ndim=len(dims)
    if ndim==1 and (name not in bidim_cols):
        names=get_band_arr1_colnames(name,dims,bands)
    elif ndim==1 and (name in bidim_cols):
        names=get_band_arr2_colnames(name,[np.sqrt(dims),np.sqrt(dims)],bands)
    elif ndim==2:
        names=get_band_arr2_colnames(name,dims,bands)
    else:
        raise ValueError("only support 1 and 2 d arrays")

    return names

def get_band_arr1_colnames(name, dims, bands):
    """
    Get db names for an array, naming 
        name_{num}
    """
    names=[]
    for i in xrange(dims[0]):
        n=bands[i]
        names.append( '%s_%s' % (name,n) )

    return names

def get_band_arr2_colnames(name, dims, bands):
    """
    Get db names for an array, naming 
        name_{num1}_{num2}
    """
    names=[]
    for i1 in xrange(dims[0]):
        n1=bands[i1]
        for i2 in xrange(dims[1]):
            n2=bands[i2]
            names.append( '%s_%s_%s' % (name,n1,n2) )

    return names

def get_oracle_type(nt):
    if 'f4' in nt:
        ot='binary_float'
    elif 'f8' in nt:
        ot='binary_double'
    elif 'i1' in nt or 'u1' in nt:
        ot='number(3)'
    elif 'i2' in nt or 'u2' in nt:
        ot='number(5)'
    elif 'i4' in nt:
        ot='number(10)'
    elif 'i8' in nt:
        ot='number(19)'
    elif 'u8' in nt:
        ot='number(20)'
    elif 'S' in nt:
        slen=nt[2:]
        ot='varchar2(%s)' % slen
    else:
        raise ValueError("unsupported numpy type: '%s'" % nt)

    return ot

def get_fits_type(name):
    if name == "id": 
        format = 'K'
    elif name == "number":
        format = 'J'
    elif 'nimage_tot' in name:
        format = 'J'
    elif name == "fofid":
        format = 'K'
    elif name == "ra":
        format = 'D'
    elif name == "dec":
        format = 'D'
    #elif name == "coadd_run":
    #    format = '27A'
    elif name == "flags":
        format = 'J'
    elif name == "time_last_fit":
        format = 'D'
    elif name == "box_size":
        format = 'I'
    elif name == "obj_flags": 
        format = 'J'
    elif "psf_flags" in name: 
        format = 'J'
    elif "psf_flux" in name:
        format = 'D'
    elif "psf_flux_err" in name: 
        format = 'D'
    elif "psf_flux_s2n" in name: 
        format = 'D'
    elif "psf_mag" in name: 
        format = 'D'
    elif "nimage_use" in name: 
        format = 'J'
    elif name == "mask_frac": 
        format = 'D'
    elif name == "psfrec_T": 
        format = 'D'
    elif "psfrec_g" in name: 
        format = 'D'
    elif name == "cm_flags": 
        format = 'J'
    elif "cm_pars" in name: 
        format = 'D'
    elif "cm_pars_cov" in name: 
        format = 'D'
    elif name == "cm_T": 
        format = 'D'
    elif name == "cm_T_err": 
        format = 'D'
    elif name == "cm_T_s2n": 
        format = 'D'
    elif "cm_flux_cov" in name: 
        format = 'D'
    elif "cm_flux" in name: 
        format = 'D'
    elif "cm_flux_s2n" in name: 
        format = 'D'
    elif "cm_mag" in name: 
        format = 'D'
    elif "cm_logsb" in name: 
        format = 'D'
    elif "cm_g" in name: 
        format = 'D'
    elif "cm_g_cov" in name: 
        format = 'D'
    elif name == "cm_weight": 
        format = 'D'
    elif name == "cm_max_flags": 
        format = 'J'
    elif "cm_max_pars" in name: 
        format = "D"
    elif "cm_max_pars_cov" in name: 
        format = 'D'
    elif name == "cm_max_T": 
        format = 'D'
    elif name == "cm_max_T_err": 
        format = 'D'
    elif name == "cm_max_T_s2n": 
        format = 'D'
    elif "cm_max_flux_cov" in name: 
        format = "D"
    elif "cm_max_flux" in name: 
        format = 'D'
    elif "cm_max_flux_s2n" in name: 
        format = 'D'
    elif "cm_max_mag" in name: 
        format = 'D'
    elif "cm_max_logsb" in name: 
        format = 'D'
    elif name == "cm_s2n_w": 
        format = 'D'
    elif name == "cm_chi2per": 
        format = 'D'
    elif name == "cm_dof": 
        format = 'D'
    elif name == "cm_flags_r": 
        format = 'J'
    elif name == "cm_s2n_r": 
        format = 'D'
    elif name == "cm_T_r": 
        format = 'D'
    elif name == "cm_psf_T_r": 
        format = 'D'
    elif name == "cm_fracdev": 
        format = 'E'
    elif name == "cm_fracdev_noclip": 
        format = 'E'
    elif name == "cm_fracdev_err": 
        format = 'E'
    elif name == "cm_TdByTe": 
        format = 'E'
    elif name == "cm_TdByTe_noclip": 
        format = 'E'
    elif name == "cm_mof_flags": 
        format = 'J'
    elif name == "cm_mof_num_itr": 
        format = 'J'
    elif "cm_mof_abs_diff" in name: 
        format = 'D'
    elif "cm_mof_frac_diff" in name: 
        format = 'D'
    elif "cm_mof_err_diff" in name: 
        format = 'D'
    elif name == 'fofind': 
        format = 'K'

    return format

def check_name(name,cols,prev,strip,printout=False):
    check = False
    n = 0
    colname = prev
    for enum,col in enumerate(cols):
        if printout:
            print col,name[0:len(name)-strip]
        if col == name[0:len(name)-strip]:
            check = True
            colname = col
            break        
    return (check,colname)

def build_table(mof_dir,mof_file,coaddid_dir):
    mof_hdu = fits.open(mof_dir+mof_file)
    mof_tab = mof_hdu[1].data
    coaddid_hdu = fits.open(os.path.join(coaddid_dir,mof_file[:12]+'_ids.fits'))
    coaddid_tab = coaddid_hdu[1].data

    defs = {}
    descr = mof_tab.view(np.ndarray).dtype.descr
    alldefs = get_coldefs(descr,defs,bands,band_cols,bidim_band_cols) 
    names = [d[0] for d in alldefs]
    formats = [d[2] for d in alldefs]
    cols = []
    prev_colname = ''
    prev_colname_band = ''
    prev_colname_multi = ''
    prev_colname_bidim = ''
    prev_colname_bidim_band = ''
    k = 0
    m = 0
    for i in xrange(len(names)):
	strip = len(names[i].split('_')[len(names[i].split('_'))-1])
        check_cols,colname = check_name(names[i],onedim_cols,prev_colname,0)
        check_band_cols,colname_band = check_name(names[i],band_cols,prev_colname_band,strip+1)
        check_bidim_cols,colname_bidim = check_name(names[i],bidim_cols,prev_colname_bidim,strip+1)
        check_bidim_band_cols,colname_bidim_band = check_name(names[i],bidim_band_cols,prev_colname_bidim_band,4)
        check_multi_cols,colname_multi = check_name(names[i],multi_cols,prev_colname_multi,strip+1)
        if i == 0:
            n = 0
            m = 0
        if i > 0 and (prev_colname != colname or colname_band != prev_colname_band or colname_bidim != prev_colname_bidim or colname_bidim_band != prev_colname_bidim_band or colname_multi != prev_colname_multi):
            n = 0
            m = 0
            k = k+1
        if check_band_cols == True:
            cols.append(fits.Column(name=names[i],format=get_fits_type(names[i]),array=mof_tab[colname_band][:,n]))
            n = n + 1
            prev_colname_band = colname_band
        elif check_bidim_cols == True:
            cols.append(fits.Column(name=names[i],format=get_fits_type(names[i]),array=mof_tab[colname_bidim][:,n,m]))
            if n == len(mof_tab[colname_bidim][0])-1:
                n = 0
                m = m + 1
            else:
                n = n + 1
            prev_colname_bidim = colname_bidim
        elif check_bidim_band_cols == True:
            cols.append(fits.Column(name=names[i],format=get_fits_type(names[i]),array=mof_tab[colname_bidim_band][:,n,m]))
            if n == len(mof_tab[colname_bidim_band][0])-1:
                n = 0
                m = m + 1
            else:
                n = n + 1
            prev_colname_bidim_band = colname_bidim_band
        elif check_multi_cols == True:
            cols.append(fits.Column(name=names[i],format=get_fits_type(names[i]),array=mof_tab[colname_multi][:,n]))
            if n == len(mof_tab[colname_multi][0])-1:
                n = 0
            else:
                n = n + 1
            prev_colname_multi = colname_multi
        else:
            if names[i] == "id":
                newname = "coadd_object_id"
            elif names[i] == "number":
                newname = "nbr"
            else:
                newname = names[i]
            cols.append(fits.Column(name=newname,format=get_fits_type(names[i]),array=mof_tab[names[i]]))
            prev_colname = colname
    #print(len(mof_tab))
    tilename_name = [mof_file[:12]] * len(mof_tab)	
    cols.append(fits.Column(name='TILENAME',format='27A',array=tilename_name))
    for b,band in enumerate(bands):
	psf_mag_err = 1.0857362*mof_tab['psf_flux_err'][:,b]/mof_tab['psf_flux'][:,b]
   	psf_mag_err[psf_mag_err == 1.0857362] = -9999
        cm_mag_err = 1.0857362*np.sqrt(mof_tab['cm_flux_cov'][:,b,b])/mof_tab['cm_flux'][:,b]
        cm_mag_err[cm_mag_err == 1.0857362] = -9999
        cm_mag_err = [-9999 if np.isnan(x) else x for x in cm_mag_err]
    	cols.append(fits.Column(name='psf_mag_err_'+band,format='D',array=psf_mag_err))
        cols.append(fits.Column(name='cm_mag_err_'+band,format='D',array=cm_mag_err))
    new_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    mask = np.in1d(new_hdu.data['coadd_object_id'],coaddid_tab['COADD_OBJECT_ID'])
    print 'Before masking:',len(new_hdu.data),'After masking:',len(new_hdu.data[mask])
    new_tbdata = new_hdu.data[mask]
    new_hdu = fits.BinTableHDU(data=new_tbdata)

    return new_hdu

def main():

    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--datadir",dest="mof_dir",help="Directory where data is read from",default='/home/nsevilla/y3mof/test/tilelinks/')
    parser.add_option("-i","--iddir",dest="coaddid_dir",help="Directory where coadd object ids are stored per tile",default='/home/nsevilla/y3mof/test/mof_ids_test/')
    parser.add_option("-o","--outdir",dest="out_dir",help="Directory where data for upload will be stored",default='/home/nsevilla/y3mof/test/mof_uploads/')

    # Parse command line
    (options, args) = parser.parse_args()

    if not os.path.isdir(options.mof_dir):
        print 'Path with data not found at',options.mof_dir
        sys.exit(1)
    if not os.path.isdir(options.coaddid_dir):
        print 'Path with coadd id files not found at',options.coaddid_dir
        sys.exit(1)
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    print 'Reading data from',options.mof_dir
    print 'Reading ids from',options.coaddid_dir
    print 'Writing files to upload to',options.out_dir

    #get list of files
    print 'Getting list of files...'
    mof_files = [f for f in os.listdir(options.mof_dir) if os.path.isfile(os.path.join(options.mof_dir, f))]

    #mof_files =  ['DES0051+0252-y3v02-mof-001.fits']

    for f,mof_file in enumerate(mof_files):
        if "DES" in mof_file: 
            print 'Processing tile',f+1,'/',len(mof_files),mof_file
            if not os.path.isfile(options.mof_dir+mof_file):
                print 'File',options.mof_dir+mof_file,'not found'
                sys.exit()
            out_hdu = build_table(options.mof_dir,mof_file,options.coaddid_dir)
            #if f == 0:
            #    out_hdu = new_hdu
            #else:
            #    for colname in out_hdu.data.columns.names:
            #        out_hdu.data[colname][out_hdu.data.shape[0]:] = new_hdu.data[colname]
        out_file = os.path.join(options.out_dir,mof_file[0:12]+'_mof_upload.fits')
        out_hdu.writeto(out_file,clobber='True') 

if __name__ == '__main__':
    sys.exit(main())


