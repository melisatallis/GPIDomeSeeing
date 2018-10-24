import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.fftpack as fft
from astropy.io import fits
from poppy import zernike
import gpi_phase_analysis as gpa
import sys
import os

def main(filename, save_images=True):
    """
    Melisa Tallis - 2018-10-10
    Do a spatial frequency power analysis of GPI phase and save data products to specific directory.
    """
    
    ## Gemini telescope parameters
    outD = 7.77010  # primary diameter (m)
    inD = 1.024     # inner M2 diameter (m)
    
    ## GPI DM parameters
    n = 48.           # number sample points across the screen 
    nacross = 43.2    # number of subapertures across the aperture
    
    ## phase sample parameters
    pscale = outD/(nacross)     #  pixel size (m) of samples in pupil plane

    ## make the aperture
    ap = gpa.MakeAperture(48,pscale)
    ap_nan = np.copy(ap.astype(np.float))
    ap_nan[np.where(ap == 0)] = np.nan
        
    ## open the file 
    #rootdir = '/home/sda/mtallis/GPIDomeSeeing/PhaseScripts/aotelem/Reduced/20160227/'
    rootdir = "/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/Reduced/20160227/"
    hdulist = fits.open(rootdir+filename,memmap=True)
    phase = hdulist[0].data
    
    ## Get the phase dimensions
    phdim = phase.shape # output is in 
    phx   = phdim[1]
    phy   = phdim[2]
    timesteps = phdim[0]
    
    avg_phase = np.nanmean(phase*ap_nan[None,:,:],axis=0)
    
    ## Remove zernikes
    z_basis = zernike.zernike_basis_faster(nterms= 6, npix = 48)
    z_coeff = zernike.opd_expand_nonorthonormal(avg_phase, aperture=aperture, nterms=6)
    thin_lens = np.sum(z_coeff[:,None,None]*z_basis[:,:,:],axis=0)
    print('Done removing zernikes')

    c_phase = (phase - thin_lens)*ap_nan
    c_phase[np.isnan(c_phase)]=0.
    
    ## Compute the DFT of phase
    phFT = np.zeros((timesteps,phx,phy), dtype=complex)
    for t in np.arange(timesteps):
        phFT[t,:,:] = fft.fftshift(fft.fft2(c_phase[t,:,:]*ap))/(ap.sum)
    print('Done with FT')
    
    ## Compute the PSD of phase
    psd2D = np.zeros((timesteps, phx, phy),dtype=float)
    for k in np.arange(phx):
        for l in np.arange(phy):
            psd2D[:,k,l] = np.abs(phFT[:,k,l])**2
    
    print('Done with PSD')
    avg_psd2D = np.mean(psd2D, axis=0)
    
    ## Create k-grid
    kr = gpa.MakeFreqGrid(phx,pscale) 
    
    ## Compute radial slices of 2D PSD
    avg_psd1D =  gpa.RadialProfile(avg_psd2D)
    freq = gpa.RadialProfile(kr)
    
    # Fit a power law to the 1D PSD
    slope = gpa.PowerLawFit(avg_psd1D,freq)
    print(slope)   
    if save_images:
        gpa.PSDplot(avg_psd1D,freq)
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S')
        psdfileroot = (filename + 'sp_psd'+ timestamp + '.png')   #  Need to parse out date from filename
        outdir = "/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/sp_psd"
        psdoutfile = os.path.join(outdir, psdfileroot)
        plt.savefig(psdoutfile)
          
"""if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help='Name of data file.')
    parser.add_argument("directory", help='Main directory to perform all analysis')

    args = parser.parse_args()
    directory = args.directory
    filename = args.filename
    
    main(filename, outdir)        
    """
     
        
        
    
    
    