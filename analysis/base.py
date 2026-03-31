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

flag_info = [
    (1, "DO_NOT_USE", 0),
    (2, "SATURATED", 1),
    (4, "JUMP_DET", 2),
    (8, "DROPOUT", 3),
    (16, "OUTLIER", 4),
]

def check_flagged_pix(all_dq):
    ''' check what pixels are flagged in a jwst fits file. '''
    
    unique_dq = np.unique(all_dq)
    print(f"\nunique GROUPDQ values in data: {len(unique_dq)}")
    if len(unique_dq) <= 30:
        print(f"all unique values: {unique_dq}")
    else:
        print(f"sample values: {unique_dq[:30]}")
    dq_flat = all_dq.flatten()
    n_total = len(dq_flat)
    print(f"\nGROUPDQ flag statistics:")
    print(f"  DQ = 0 (good): {np.sum(dq_flat == 0)} ({100*np.sum(dq_flat == 0)/n_total:.2f}%)")
    print("\nflag breakdown:")
    for flag_val, flag_name, bit_num in flag_info:
        n_flagged = np.sum((dq_flat & flag_val) > 0)
        if n_flagged > 0:
            print(f"  Bit {bit_num} - {flag_name} (value {flag_val}): {n_flagged} ({100*n_flagged/n_total:.3f}%)")

def mask_bad_pixels(all_dq):
    ''' create a good pixel mask using jwst's groupdq flags. '''
    
    good_pixel_mask = (all_dq == 0)
    n_good = np.sum(good_pixel_mask)
    n_total = good_pixel_mask.size
    print(f"\nmask using GROUPDQ flags (GROUPDQ == 0):")
    print(f"  good pixels: {n_good} ({100*n_good/n_total:.2f}%)")
    
    return good_pixel_mask

def make_global_background(frame, box_size=64, filter_size=3): 
    ''' simple global background estimator for fits image data. '''
    
    sigma_clip = SigmaClip(sigma=3.0)                                    
    bkg = Background2D(frame, box_size=box_size,                         
                        filter_size=filter_size,                          
                        sigma_clip=sigma_clip,                            
                        bkg_estimator=MedianBackground())      
    
    return bkg.background

def aperture_photometry_zeroframe(frame, x_cent, y_cent,                                                                            
                                    tframe,                                                                                           
                                    ap_radius=2.5,                                                                                    
                                    ann_inner=4.5,                                                                                    
                                    ann_outer=10.5,
                                    sigma=2.5,                                                                                        
                                    maxiters=5):
      """                                                                                                                             
      aperture photometry for zeroframe data already in e-/s.

      Parameters
      ----------
      frame    : 2D array in e-/s (gain, tframe, and flat already applied)
      tframe   : float, frame time in seconds (for photon noise calculation)                                                          
      """                                                
                                        
      y_grid, x_grid = np.indices(frame.shape)                                                                                        
      r = np.sqrt((x_grid - x_cent)**2 + (y_grid - y_cent)**2)                                                                        
                                                                                                                                      
      aperture_mask = (r <= ap_radius) & np.isfinite(frame)                                                                           
      annulus_mask  = (r >= ann_inner) & (r <= ann_outer) & np.isfinite(frame)                                                        
                                                                                                                                      
      # background in e-/s per pixel
      _, median_bkg, std_bkg = sigma_clipped_stats(                                                                                   
          frame[annulus_mask], sigma=sigma, maxiters=maxiters                                                                         
      )                                                                                                                               
                                                                                                                                      
      # net flux in e-/s                                                                                                              
      n_aperture_pixels = np.sum(aperture_mask)                                                                                       
      aperture_sum = np.sum(frame[aperture_mask])
      ap_area = np.pi * ap_radius**2
      flux_rate = aperture_sum - median_bkg * ap_area                                                                       
                                                                                                                                      
      # noise: compute in electron space, convert back to rate                                                                        
      flux_electrons = flux_rate * tframe                                                                                             
      std_bkg_e = std_bkg * tframe                                                                                                    
                  
      n_clipped = np.sum(
          np.abs(frame[annulus_mask] - median_bkg) <= sigma * std_bkg
      )                                                                                                                               
                                                                                                                                      
      variance_electrons = (                                                                                                          
          np.abs(flux_electrons) +                    # Poisson noise                                                      
          ap_area * std_bkg_e**2 +                    # background noise
          ap_area**2 * std_bkg_e**2 / n_clipped       # background uncertainty                                              
      )
      flux_error_rate = np.sqrt(variance_electrons) / tframe                                                                          
                                                                                                                                      
      return {                                                                                                                        
          'count_rate':       flux_rate,            # e-/s                                                                            
          'count_rate_error': flux_error_rate,      # e-/s
          'background_rate':  median_bkg,           # e-/s per pixel                                                                  
          'exposure_time':    tframe,
          'n_aperture_pixels': n_aperture_pixels,                                                                                     
          'n_clipped_pixels':  int(n_clipped)                                                                                         
      }     

def zeroframe_corr_pipeline(cleaned_ramp_filenames, calibrated_filenames, dir):
    ''' convert zeroframes to e-/s. '''
    
    cleaned_zeroframes = extract_zeroframes_from_ramp(cleaned_ramp_filenames)
    
    with fits.open(cleaned_ramp_filenames[0]) as hdul:
        zhdr_c = hdul['ZEROFRAME'].header
        shdr_c = hdul['SCI'].header
        ghdr_c = hdul[0].header

    with fits.open(calibrated_filenames[0]) as hdul:
        shdr_cal = hdul['SCI'].header
        ghdr_cal = hdul[0].header

    zeroframe_exp = ghdr_c['TFRAME']

    # get nrcb4 noise information
    gain_nrcb4_fn = get_ref_file(ghdr_c,   'R_GAIN', dir) 
    flat_nrcb4_fn = get_ref_file(ghdr_cal, 'R_FLAT', dir)

    with fits.open(flat_nrcb4_fn) as hdul:                                                                                              
        flat_data = hdul['SCI'].data.copy()                                                                        
        flat_dq   = hdul['DQ'].data.copy()         

    with fits.open(gain_nrcb4_fn) as hdul:        
        gain_data = hdul['SCI'].data.copy()

    # convert zeroframe to e-/s using the gain file + flat field the data
    good_gain = (gain_data > 0) & np.isfinite(gain_data)                                                                                
    good_flat = (flat_dq == 0) & (flat_data > 0.1)
    good = good_gain & good_flat                                                                                                        

    zf_rate = np.full_like(cleaned_zeroframes, np.nan, dtype=np.float64)                                                    
    for i in range(cleaned_zeroframes.shape[0]):
        for j in range(cleaned_zeroframes.shape[1]):                                                                              
            frame = cleaned_zeroframes[i, j].astype(np.float64)
            frame[good] /= gain_data[good]       # DN → electrons                                                                       
            frame[good] /= zeroframe_exp         # electrons → e-/s                                                                     
            frame[good] /= flat_data[good]       # flat field                                                                           
            frame[~good] = np.nan                                                                                                       
            zf_rate[i, j] = frame

    return zf_rate
    
