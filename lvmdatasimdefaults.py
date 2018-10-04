param = {
    "verbose":"yes",
    "minWave":3550.0,
    "maxWave":9850.0,
    "waveStep":0.1,
    "sky_conditions":"dark",
    "seeing_fwhm_ref":1.1,
    "seeing_wlen_ref":6355,
    "seeing_moffat_beta":3.5,
    "moon_phase":0.5,
    "moon_zenith":100,
    "separation_angle":60,
    "airmass":1.0,
    "telescope_name":"LVM",
    "tel_mirror_diameter":1.0,
    "tel_obscuration_diameter":0.47,
    "tel_support_width":0.08,
    "tel_fiber_diameter":107.0,
    "tel_field_radius":414.0,
    "fiberloss_method":"perfect",
    "instr_sigma1d":5.1,
    "b_read_noise":3.0, # electron/pixel**2
    "b_dark_current":3.0, # electron/(hour pixel**2)
    "b_gain":1.0, # electron/adu
    "b_output_pixel_size":0.5,# Angstrom
    "b_psf":"specpsf/psf-quicksim.fits",
    "b_throughput":"throughput/thru-b.fits",
    "r_read_noise":2.9, # electron/pixel**2
    "r_dark_current":2.0, # electron/(hour pixel**2)
    "r_gain":1.0, # electron/adu
    "r_psf":"specpsf/psf-quicksim.fits",
    "r_throughput":"throughput/thru-r.fits",
    "z_read_noise":2.9, # electron/pixel**2
    "z_dark_current":2.0, # electron/(hour pixel**2)
    "z_gain":1.0, # electron/adu
    "z_output_pixel_size":0.5, # Angstrom
    "z_psf":"specpsf/psf-quicksim.fits",
    "z_throughput":"throughput/thru-z.fits",
    "observatory":"APO",
    "exposure_time":1000.0, # s
    "temperature":15, # deg_C
    "realative_humidity":0., #
    "exposure_start_time":55000.5
}
