import numpy as np
import scipy.fftpack as fft

def makeFreqGrid(n,pscale):
    kx = fft.fftshift(fft.fftfreq(n,pscale))
    ky = fft.fftshift(fft.fftfreq(n,pscale))
    mg = np.meshgrid(kx,ky)
    kr = np.sqrt(np.sum((m**2 for m in mg))) 
    return kr

def radialProfile(image, center=None):
    """
    Calculate the avearge radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    ## Calculate the indices from the image
    y,x = np.indices((image.shape)) # first determine radii of all pixels
    
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
     
    r = np.hypot(x - center[0], y - center[1]).astype(np.int) 
    
    n = np.bincount(r.ravel())
    sy = np.bincount(r.ravel(), image.ravel())
    sy2 = np.bincount(r.ravel(), image.ravel()*image.ravel())
    mean = sy/n
    std = np.sqrt(sy2/n - mean*mean)
    
    return mean,std,n 

'''def PowerLawFit(psd,freq,low_bound = 0., up_bound = 2.7):
    """Plot the 1D psd"""
    
    ##  create matplotlib figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    def func(x, a, b):
        return a+(b*x) 
    
    par,pcov = optimize.curve_fit(func,np.log10(freq[(freq > low_bound) & (freq < up_bound)]),
                                  np.log10(psd[(freq > low_bound) & (freq < up_bound)]), p0=(1, -3.5))
    slope = par[0]
    intercept = par[1]
    return slope

def PSDplot(psd, freq, low_bound = 0., up_bound = 2.7, display=False, **kwargs):
    """Plot the 1D psd"""
    
    ##  create matplotlib figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    def func(x, a, b):
        return a+(b*x) 
    
    par,pcov = optimize.curve_fit(func,np.log10(freq[(freq > low_bound) & (freq < up_bound)]),
                                  np.log10(psd[(freq > low_bound) & (freq < up_bound)]), p0=(1, -3.5))
    slope = par[0]
    intercept = par[1]
    
    ## Plot original PSD and linear fit
    img = ax.loglog((freq[(freq > low_bound) & (freq < up_bound)]),(psd[(freq > low_bound) & (freq < up_bound)]),'b.',(                     freq[(freq > low_bound) & (freq < up_bound)]),
                    10**(func(np.log10(freq[(freq > low_bound) & (freq < up_bound)]),*par)), 'r')

    ax.legend(['PSD', 'slope = {0:.2f}, intercept={1:.2f}'.format(slope, intercept)],loc=3, fontsize=10)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='grey', linestyle='-')
    ax.set_ylabel('Power Spectrum')
    ax.set_xlabel('Spatial Frequency [1/m]')
    
    ## Modify plot based on keyword arguments
    if 'title' in kwargs: ax.set_title(kwargs['title'], fontsize=10)
    if 'xlabel' in kwargs: ax.set_xlabel(kwargs['xlabel'], fontsize=15)
    if 'ylabel' in kwargs: ax.set_ylabel(kwargs['ylabel'], fontsize=15)
    if display: plt.show()
    if 'save_image' in kwargs: plt.savefig(kwargs['save_image'])    
    
    return slope'''
