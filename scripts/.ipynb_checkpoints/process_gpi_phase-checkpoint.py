import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.fftpack as fft
from astropy.io import fits
from poppy import zernike
import gpi_phase_analysis as gpa
import sys
import os

#  Physical Parameters
outD = 7.77010            # primary diameter (m)
inD = 1.024               # inner M2 diameter (m)
n = 48                    # number sample points across the screen (Not the number of subapertures)
nacross = 43              # number of subapertures across the aperture
pscale = outD/(nacross)   # pixel size (m) of samples in pupil plane
     
        
#  Make the aperture 
ap = gpa.makeAperture(n,pscale)         #contains zeros
ap_nan = np.copy(ap.astype(np.float))  
ap_nan[np.where(ap==0)] = np.nan        # contains nans       
    
#  Make frequency grid    
kr = gpa.makeFreqGrid(n,pscale)    

#  Search directory tree for dm phase cubes
#rootdir = "/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/Reduced/"
fname_list = list()
name_list = list()

for root, dirs, files in os.walk(rootdir):
    for name in files:
        (base,ext) = os.path.splitext(name)
        if ext in ('.fits') and base[0]=='a':
            full_name = os.path.join(root,name)
            fname_list.append(full_name)  
            name_list.append(name)

#  Store filename and outputs inside dataframe
df = pd.DataFrame({'path':[],'filename':[],'whenstr':[],'slope':[]})
df['path'] = fname_list
df['filename'] = name_list

#  begin psd analysis
n=0
for file in df.loc[0:5,'path']:
    
    hdulist = fits.open(file,memmap=True)        # phase data shape (time,xpix,ypix)
    df.loc[n,'whenstr'] = hdulist[0].header['whenstr']  
    
    phase = hdulist[0].data.astype('float')
    timesteps, phx, phy = phase.shape            # contains a datacube
    avg_phase = np.nanmean(phase*ap_nan,axis=0)  # used to find average zernikes 
    
    # remove zernikes form cube
    z_basis = zernike.zernike_basis_faster(nterms= 6, npix = 48)
    z_coeff = zernike.opd_expand_nonorthonormal(avg_phase, aperture=ap, nterms=6)
    thin_lens = np.sum(z_coeff[:,None,None]*z_basis[:,:,:],axis=0)

    c_phase = (phase - thin_lens[None,:,:])*ap_nan
    c_phase[np.isnan(c_phase)]=0.
    print('Done removing zernikes')
    
    # computer fourier transform of cube
    phFT = np.zeros((timesteps,phx,phy), dtype=complex)
    for t in np.arange(timesteps):
        phFT[t,:,:] = fft.fftshift(fft.fft2(c_phase[t,:,:]))/ap.sum()
    print('Done with FT')
    
    
    # compute 2d psd cube
    psd2D = np.zeros((timesteps, phx, phy),dtype=float)
    for k in np.arange(phx):
        for l in np.arange(phy):
            psd2D[:,k,l] = np.abs(phFT[:,k,l])**2
    avg_psd2D = np.mean(psd2D, axis=0)
    print('Done with PSD')    
    
    # compute radial average of 2d psd cube and frequency
    avg_psd1D =  radialProfile(avg_psd2D)
    freq = radialProfile(kr)
    
    outdir = '/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/sp_psd/'
    #outdir = '/Users/melisatallis/Documents/Research/GPIDomeSeeing/data/aotelem/sp_psd/'
    date = df.loc[n,'whenstr']
    # Plot spectrum
    #slope, intercept = PSDplot(avg_psd1D, freq, low_bound = 0.3, display=False, whenstr=date, filebase=outdir)
    
    low_bound = 0.
    up_bound = 3.
    whenstr = df.loc[n,'whenstr']
    
    # Plotting PSD
    fig = plt.figure(figsize=[10,7])
    ax = fig.add_subplot(111)

    def func(x, a, b):
        return a+(b*x) 
    
    par,pcov = optimize.curve_fit(func,np.log10(freq[(freq > low_bound) & (freq < up_bound)]),
                                  np.log10(avg_psd1D[(freq > low_bound) & (freq < up_bound)]), p0=(1, -3.5))
    slope = par[0]
    intercept = par[1]
    
    ## Plot original PSD and linear fit
    img = ax.loglog((freq[(freq > low_bound) & (freq < up_bound)]),(avg_psd1D[(freq > low_bound) & (freq < up_bound)]),
                    'bo',(freq[(freq > low_bound) & (freq < up_bound)]),
                    10**(func(np.log10(freq[(freq > low_bound) & (freq < up_bound)]),*par)), 'r')

    ax.legend(['PSD', 'slope = {0:.2f}, intercept={1:.2f}'.format(slope, intercept)],loc=3, fontsize=15)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='black', linestyle='-') 
    ax.set_title(whenstr, fontsize=30, y=1.04)
    ax.set_ylabel('Power Spectrum',fontsize=15)
    ax.set_xlabel('Spatial Frequency',fontsize=15)
    
    filename = "{0}_PSD{1}.png".format(outdir, whenstr)
    plt.savefig(filename)
    #plt.close(fig)
    
    df.loc[n,'slope'] = slope
    n=n+1