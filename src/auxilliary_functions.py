# Auxilliary functions for rerunning DALEC locally following calibration with CARDAMOM
import numpy as np
from scipy.stats import gaussian_kde

def pdf(serie,nbins=100,zeros=True):
    'returns pdf (x,y) from a serie with nbins'
    k = gaussian_kde(serie) #kernel density
    m = k.dataset.min()
    M = k.dataset.max()
    x = np.arange(m,M,(M-m)/nbins)
    v = k.evaluate(x) #density curve
    if zeros:
        inter=x[1]-x[0]
        x=np.array([x[0]-inter]+list(x)+[x[-1]+inter])
        v=np.array([0]+list(v)+[0])
    return {'x':x,'v':v}

def GR(chains):
    
    n = len(chains[0]) # number of draws per chain
    m = len(chains)    # number of chains
 
    # calculate the within chain variance W: mean of the variance of each chain
    W = 0.
    for cc in range(m):

        sj2 = np.var(chains[cc],ddof=1)
        W += sj2
    W = W / m

    # calculate the between chain variance B: variance of chains means
    theta = 0.
    for cc in range(m):
        theta += np.mean(chains[cc])
    theta = theta*1./m

    B = 0.
    for cc in range(m):
        B += (np.mean(chains[cc])-theta)**2
    B = B*n/(m+1)

    # estimate variance
    vartheta = (1.-1./n)*W+B/n

    return np.sqrt(vartheta/W)
