import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
from astropy.io import fits

def main(filename, outdir, save_images, aperture=True, blackman=True):
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
    N = 48.          # number sample points across the screen 
                    # (Not the number of subapertures across the aperture which is less) 
    nacross = 43.2    # number of subapertures across the aperture
    
    ## phase sample parameters
    pscale = outD/(nacross)     #  pixel size (m) of samples in pupil plane

    if aperture:
        x         = np.linspace(-(N-1)/2,(N-1)/2,N)*pscale 
        y         = np.linspace(-(N-1)/2,(N-1)/2,N)*pscale
        ax,ay     = np.meshgrid(x,y)
        ar        = np.sqrt(ax**2 +ay**2)  
        ap_outer  = (ar <= outD/2)
        ap_inner  = (ar <= inDs/2)   
        a        = (ap_outer ^ ap_inner).astype(int)
        
    # open the file 
    hdulist = fits.open(rootdir+filename,memmap=True)
    data = hdulist[0].data
    
    # Get the phase dimensions
    phdim = hdulist[0].data.shape # output is in 
    phx   = phdim[1]
    phy   = phdim[2]
    timesteps = phdim[0]
    
    phFT = np.zeros((timesteps,phx,phy), dtype=complex)
    for t in np.arange(timesteps):
        phFT[t,:,:] = fft.fftshift(fft.fft2(data[t,:,:]*a))/(a.sum)
    print('Done with FT')
    
    print 'Doing PSD'
    psd2D = np.zeros((timesteps, phx, phy),dtype=float)
    for k in np.arange(phx):
        for l in np.arange(phy):
            psd2D[:,k,l] = np.abs(phFT[:,k,l])**2
    
    varpsd = np.sum(psd, axis=0)
    
    # Create k-grid
    kx = fft.fftshift(fft.fftfreq(phx,pscale))
    ky = fft.fftshift(fft.fftfreq(phx,psacle))
    mg = np.meshgrid(kx,ky)
    kr = np.sqrt(np.sum((m**2 for m in mg))) 
    
    plt.figure(1)
    plt.loglog(kr,varpsd)
    plt.grid(True)
    plt.show()
    
     
        
        
    
    
    