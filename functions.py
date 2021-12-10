import numpy as np
from scipy import stats
import scipy
from scipy.signal import fftconvolve
from io import StringIO
from scipy.stats import gaussian_kde
from kde_diffusion import kde1d
from Cfunctions import solve_gaussian,solve_gaussianC, substractArray,sumArrays, divideArray,PDP_C,KDEp_C

#Load Data
#def loadData(filename):
#    data=np.genfromtxt(filename,skip_header=2)
#    return data
    


import pandas as pd
def loadData(filename):
    data=pd.read_csv(filename,skiprows=2,header=None,delimiter="\t",names=["Age", "Error"])
    return data


#PDP
#def solve_gaussian(x,sigma):
#    return (1. / (sigma*np.sqrt(2.*np.pi)) * np.exp(- (x) * (x) / (2. * sigma * sigma)))

def PDP(x1_grid,data_array,sigma_array):
    PDP=PDP_C(x1_grid,data_array,sigma_array)
    return PDP

#KDE- Fixed -Method 1
#Implemented from the equations presented in:
#Vermeesch, P. (2012). On the visualisation of detrital age distributions.
#Chemical Geology, 312, 190-194.

def KDEp(x1_grid,data_array,bandwidth):
    KDE=KDEp_C(x1_grid,data_array,bandwidth)
    return KDE

#def KDEp(x1_grid,data_array,bandwidth):
#
##
##
##	for i in range(0,len(data_array)):
##        x=substractArray(x1_grid,data_array[i])
##        if i==0:
##			KDE=solve_gaussian(x,bandwidth)
##		else:
##			KDE+=solve_gaussian(x,bandwidth)
#
#	return 0

#KDE-Fixed- Method 2 (Not Implemented here)
def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    #Kernel Density Estimation with Scipy
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

#Python re-implementation of the original MatLab script that describes the method of
#Botev et al., 2010.
#Python Implementation performed by John Henning -
#John Hennig. (2021, April 6). John-Hennig/KDE-diffusion: KDE-diffusion 1.0.3 (Version v1.0.3).
#Zenodo. http://doi.org/10.5281/zenodo.4663430
def kde_difussion(x1_grid,data_array,Min,Max):
    n=len(x1_grid)
    (KDE, grid, bandwidth) = kde1d(data_array,n, limits=(Min,Max))
    return KDE,bandwidth

    

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

#KDE-Adaptative + PDP
def KDEadap_PDP(data, sigma,x1_grid,Min,Max):
    Data=[]
    KDEv=kde_difussion(x1_grid,data,Min,Max)[0]
    PDPv=PDP(x1_grid,data,sigma)
    Data.append(KDEv)
    Data.append(PDPv)
    return Data

