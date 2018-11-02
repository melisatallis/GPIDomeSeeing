import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.fftpack as fft
from astropy.io import fits
from scipy import optimize
import poppy
import gpi_phase_analysis as gpa
import sys
import os
import time

#  Physical Parameters
outD = 7.77010            # primary diameter (m)
inD = 1.024               # inner M2 diameter (m)
n = 48                    # number sample points across the screen (Not the number of subapertures)
nacross = 43              # number of subapertures across the aperture
pscale = outD/(nacross)   # pixel size (m) of samples in pupil plane
     
        
#  Make the aperture 
x = np.linspace(-(n-1)/2,(n-1)/2,n)*pscale 
y = np.linspace(-(n-1)/2,(n-1)/2,n)*pscale
mg = np.meshgrid(x,y)
ar = np.sqrt(np.sum((m**2 for m in mg))) 
ap_outer = (ar <= outD/2)
ap_inner = (ar <= inD/2)   
ap = (ap_outer ^ ap_inner).astype(int)        #contains zeros

ap_nan = np.copy(ap.astype(np.float))  
ap_nan[np.where(ap==0)] = np.nan        # contains nans       
    
#  Make frequency grid    
kr = gpa.makeFreqGrid(n,pscale)    

#  Search directory tree for dm phase cubes
#rootdir = "/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/Reduced/"

rootdir = "/home/sda/mtallis/PhaseScripts/aotelem/Reduced/"
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

#  Results path
#save_path = '/Users/MelisaT/Documents/Research/GPIDomeSeeing/data/'

save_path = '/home/sda/mtallis/Results/'
dstr = time.strftime('%Y%m%d')

#  begin psd analysis
n=0
for file in df.loc[0:2,'path']:
    
    hdulist = fits.open(file,memmap=True)        # phase data shape (time,xpix,ypix)
    df.loc[n,'whenstr'] = hdulist[0].header['whenstr']  
    whenstr = df.loc[n,'whenstr']
    
    phase = hdulist[0].data.astype('float')
    timesteps, phx, phy = phase.shape            # contains a datacube
    avg_phase = np.nanmean(phase*ap_nan,axis=0)  # Identify long timescale aberrations 
    
    # remove zernikes from cube
    z_basis = poppy.zernike.zernike_basis_faster(nterms= 6, npix = 48)
    z_coeff = poppy.zernike.opd_expand_nonorthonormal(avg_phase, aperture=ap, nterms=6)
    thin_lens = np.sum(z_coeff[:,None,None]*z_basis[:,:,:],axis=0)

    c_phase = (phase - thin_lens[None,:,:])*ap_nan
    c_phase[np.isnan(c_phase)]=0.
    print('Done removing zernikes')
    
    # computer normalized DFT of cube
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
    
    # compute radial average of 2d psd cube and radial frequency
    avg_psd1D =  gpa.radialProfile(avg_psd2D)
    freq = gpa.radialProfile(kr)
    
    #  plot PSD vs. freq
    low_bound = 1/(nacross*pscale)    #  set by aperture size
    up_bound = 1/(2*pscale) #  nyquist limit
    
    fig = plt.figure(figsize=[9,6])
    ax = fig.add_subplot(111)

    def func(x, a, b):
        return a+(b*x) 
    
    par,pcov = optimize.curve_fit(func,np.log10(freq[(freq > 4*low_bound) & (freq < up_bound)]),
                                np.log10(avg_psd1D[(freq > 4*low_bound) & (freq < up_bound)]), p0=(1, -3.5))
    slope = par[0]
    intercept = par[1]
    
    ## Plot original PSD and linear fit
    img = ax.loglog((freq),(avg_psd1D),'bo',(freq[(freq > 4*low_bound) & (freq < up_bound)]),
                    10**(func(np.log10(freq[(freq > 4*low_bound) & (freq < up_bound)]),*par)), 'r',lw=2)

    ax.legend(['_nolegend_', 'slope = {0:.2f}, intercept={1:.2f}'.format(slope, intercept)],loc=3, fontsize=15)
    ax.tick_params(axis='both', which='both', labelsize=16,direction='in')
    ax.grid(b=True, which='major', color='k', linestyle='-') 
    ax.set_ylabel('Power Spectrum',fontsize=20)
    ax.set_xlabel('Spatial Frequency',fontsize=20)
    ax.set_title(whenstr, fontsize=20, y=1.04)
   
    # store plots and tables in date of observation directory (be more descriptive)
    date_dir = os.path.join(save_path,whenstr[0:8])
    gpa.make_dir(date_dir)
    
    sp_psd_dir = os.path.join(date_dir,'sp_psd')
    gpa.make_dir(sp_psd_dir)
    
    plots_dir = os.path.join(sp_psd_dir,'plots')
    gpa.make_dir(plots_dir)
    
    tables_dir = os.path.join(sp_psd_dir,'tables')
    gpa.make_dir(tables_dir)
    
    table_prefix = 'sp_psd_table_'
    plots_prefix = 'sp_psd_plot_'
    
    table = pd.DataFrame({'freq':freq, 'sp_psd':avg_psd1D})
    table.to_csv(sp_psd_dir +'/tables/'+ table_prefix + whenstr +'_' + dstr + '.csv',index=False) 
 
    plt.savefig(sp_psd_dir +'/plots/'+ plots_prefix + whenstr + '_' + dstr + '.pdf')  ## make the date input/comment 
    print('done plotting')
    
    # Store measured slope in summary dataframe
    df.loc[n,'slope'] = slope
    n=n+1

summary_prefix = 'sp_psd_summary_'
df.loc[:,('whenstr','slope')].to_csv(save_path + summary_prefix + dstr+'.csv')
    
