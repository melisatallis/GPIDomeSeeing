import os
import argspace
import sys
import numpy as np
import copy

from os.path import join, splitext

def main(filename, outdir, save_images, aperture=False, blackman=True):
    """
    #################
    Melisa Tallis - 2018-10-10
    Do a spatial frequency power analysis of GPI phase and save data products to specific directory.

    Inputs:
        filename    - (str) GPI reduced data file name.
        
    Flags:
        save_images - (bool): PSD plots
        outdir      - change to save output files to a directory other than current
    """
    
    rootdir = '.' 
    
    if save_images:
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S')
        psdfileroot = (datafile + 'sp_psd'+ timestamp + )   #  Need to parse out date from filename
        psdoutfile = os.path.join(outdir, psdfileroot)
    
    ## Gemini telescope parameters
    outD = 7.77010  # primary diameter (m)
    inD = 1.024  # inner M2 diameter (m)
    
    ## GPI DM parameters
    N = 48          # number sample points across the screen 
                    # (Not the number of subapertures across the aperture which is less) 
    nacross = 43.2    # number of subapertures across the aperture
    
    ## phase sample parameters
    pscale = outD/(nacross)     #  pixel size (m) of samples in pupil plane
    T = 1.0 / pscale            #  sample spacing
    
    if aperture:
        ax, ay    = generate_grids(N, scalefac=pscale) # need to create seperate script
        ar        = np.sqrt(ax**2 + ay**2) ## aperture radius
        ap_outer  = (ar <= outD/2)
        ap_inner  = (ar <= inDs/2)   
        ap        = (ap_outer ^ ap_inner).astype(int)
        
    if blackman:
        wx, wy    = generate_grids(N, scalefac=pscale) # need to create seperate script
        wr        = np.sqrt(wx**2 + wy**2) ## window radius
        w         = np.blackman(wr)
        
    # open the file 
    hdulist = fits.open(rootdir+filename,memmap=True)
    data = hdulist[0].data
    
    # Get the phase dimensions
    phdim = hdulist[0].data.shape # output is in 
    phx   = phdim[1]
    phy   = phdim[2]
    timesteps = phdim[0]
    
    # Get the fourier transform of the phase at each timestep
    phFT = np.zeros((timesteps,phx,phy), dtype=complex)
    for t in np.arange(timesteps):
        phFT[t,:,:] = fft.fft2(hdulist[0].data[t,:,:]) / (phx*phy)
    print 'Done with FT'
    
    print 'Doing PSD'
    mft = np.sum(phFT, axis=0)
    kx, ky = gg.generate_grids(phx, scalefac=2*np.pi/(bign*pscale), freqshift=True)
        
        
    
    
    