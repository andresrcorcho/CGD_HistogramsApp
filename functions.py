import numpy as np
from scipy import stats
import scipy
from scipy.signal import fftconvolve
from io import StringIO
from scipy.stats import gaussian_kde

#Load Data
def loadData(filename):   
    data=np.genfromtxt(filename,skip_header=2)
    return data


#PDP
def solve_gaussian(x,sigma):
    return (1. / (sigma*np.sqrt(2.*np.pi)) * np.exp(- (x) * (x) / (2. * sigma * sigma)))

def PDP(x1_grid,data_array,sigma_array):
    for i in range(0,len(data_array)):
        if i==0:
            PDF=solve_gaussian((x1_grid-data_array[i]),sigma_array[i])
        else:
            PDF+=(solve_gaussian((x1_grid-data_array[i]),sigma_array[i]))
        
    return PDF/len(data_array)

#KDE- Fixed -Method 1
def KDEp(x1_grid,data_array,bandwidth):
	for i in range(0,len(data_array)):
		if i==0:
			KDE=solve_gaussian((x1_grid-data_array[i]),bandwidth)
		else:
			KDE+=solve_gaussian((x1_grid-data_array[i]),bandwidth)
	#counter=counter+1

	return KDE/len(data_array)

#KDE-Fixed- Method 2 (Implemented here)
def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    #Kernel Density Estimation with Scipy
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

#Peak detector
#Peak detector
def peakdet(v, delta, x = None):
    maxtab = []
    mintab = []
    
    if x is None:
        x = np.arange(len(v))
    
    v = np.asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    
    lookformax = True
    
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


#KDE-Fixed + PDP
def KDE_PDP(data, sigma,x_grid,bandwidth):
    Data=[]
    KDEv=KDEp(x_grid,data,bandwidth=bandwidth)
    PDPv=PDP(x_grid,data,sigma)
    Data.append(KDEv)
    Data.append(PDPv)
    return Data



