def ramp_files_per_exp(exp, dir):
    ''' get all ramp files for a given exposure. '''

    filenames = []
    for i in range(len(detector_list)):
        filename = 'jw05381003001_0210{}_00001_nrcb{}_ramp.fits'.format(exp, detector_list[i])
        filenames.append(dir + filename)
    print(f"found {len(filenames)} ramp files for exposure {exp}")

    return filenames


def ramp_files_per_det(det_no):
    ''' get all ramp files for a given detector. '''

    filenames = []
    for i in range(1, 13):
        if i >= 10:
            i = abcs[i-10]
        filename = 'jw05381003001_0210{}_00001_nrcb{}_ramp.fits'.format(i, det_no)
        filenames.append(dir + filename)
    print(f"found {len(filenames)} ramp files for detector {det_no}")

    return filenames

def load_all_ramp_data_with_dq(filenames):
    ''' get the data and dataquality flags for the given ramp files. '''

    n_files = len(filenames)
    print(f"\nloading {n_files} ramp files with dq flags")
    with fits.open(filenames[0]) as hdul:
        nint, ngroup, ny, nx = hdul["SCI"].data.shape
        print(f"\ndata shape: ({nint} int, {ngroup} groups, {ny} by {nx} pixels)") 
    all_data = np.zeros((n_files, nint, ngroup, ny, nx), dtype=np.float32)
    all_groupdq = np.zeros((n_files, nint, ngroup, ny, nx), dtype=np.uint8)
    for i, filename in enumerate(filenames):
        # if i % 3 == 0:
        #     print(f"  loading file {i+1}/{n_files}...")
        with fits.open(filename) as hdul:
            all_data[i] = hdul["SCI"].data
            all_groupdq[i] = hdul["GROUPDQ"].data
    print(f"\nloaded data shape: {all_data.shape}")
    print(f"loaded GROUPDQ shape: {all_groupdq.shape}")
    
    return all_data, all_groupdq

def extract_zeroframes_from_ramp(ramp_filenames):
    ''' extract the zeroframe data from the given ramp files. '''

    n_files = len(ramp_filenames)
    print(f"\nloading {n_files} ramp files")
    with fits.open(ramp_filenames[0]) as hdul:
        nint, ny, nx = hdul['ZEROFRAME'].data.shape
        print(f"\ndata shape: ({nint} int, {ny} by {nx} pixels)") 
    zeroframe_data = np.zeros((n_files, nint, ny, nx), dtype=np.float32)
    for i, filename in enumerate(ramp_filenames):
        with fits.open(filename) as hdul:
            zeroframe_data[i] = hdul['ZEROFRAME'].data
    print(f"\nloaded data shape: {zeroframe_data.shape}")

    return zeroframe_data

def open_fits(filename):
  ''' open a jwst fits file and return science and time information. '''
  
    with fits.open(filename) as hdul:
        sci = hdul["SCI"].data  
        sci_hdr = hdul["SCI"].header
        gen_hdr = hdul[0].header
        int_times = hdul['INT_TIMES'].data
    return sci, sci_hdr, gen_hdr, int_times

def get_ref_file(ghdr, keyword, ref_dir):
  ''' find reference file for specific keyword for jwst fits files. '''

      val = ghdr.get(keyword, '')
      filename = val.replace('crds://', '')  

      return ref_dir + filename


def plot_fits(data, label = None, title = None):
  ''' quick plotter for fits file data. '''

    plt.figure(figsize=(8, 8))
    plt.imshow(data, origin="lower", cmap="gray", vmin=np.nanpercentile(data, 1), vmax=np.nanpercentile(data, 99))
    if label:
        plt.colorbar(label=label)
    if title: 
        plt.title(title)
    plt.show()

def save_fits_cube(fits_cube, filename, reference_wcs=None):
    ''' 
    save difference cube with WCS.
    Example: 
        wcs = WCS(reference_header)
        save_fits_cube(difference_cube, 'output.fits', reference_wcs=wcs)
    '''
    cube_for_fits = np.nan_to_num(fits_cube, nan=0.0)
    print(f"\nSaving FITS cube: {filename}")
    print(f"  Shape: {cube_for_fits.shape}")
    print(f"  Size: {cube_for_fits.nbytes / 1e9:.2f} GB")
    
    hdu = fits.PrimaryHDU(cube_for_fits)

    if reference_wcs is not None:
        hdu.header.update(reference_wcs.to_header())
    
    hdu.header['NFRAMES'] = fits_cube.shape[0]
    hdu.header['BUNIT'] = 'ADU'
    
    hdu.writeto(filename, overwrite=True)
    print(f"Saved: {filename}")
