# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 10:10:54 200015

@author: christian

"""

import os, sys
from matplotlib import *
#import pickle
#import Dominance_onePop
import numpy
import re
import dadi
import itertools
import math

def readTraj(fn):
    f = open(fn, 'r')
    out = []
    for l in f:
        out.append( map(float, re.findall(r'\d+.\d+', l)) )
    f.close()
    return numpy.asarray(out)

GRIDSIZE = 200

print "Successfully loaded"
print "Current working directory: "+os.getcwd()

def parse_fold_sfs(sfs, sampleSize = 0, fold=True, maskSingletons=False):   # if fold=True, the sfs must already be a folded SFS!
    if fold==True:
        addZeroLength = sampleSize - len(sfs)
        if addZeroLength < 1: return("Wrong use of parse_fold_sfs function!")
        sfsout = [0.] + sfs + [0.]*addZeroLength
        if maskSingletons==True:
            sfsout = dadi.Spectrum(sfsout, data_folded=True, mask = [True]*2 + [False]*(len(sfs)-1) + [True]*addZeroLength)
        else:
            sfsout = dadi.Spectrum(sfsout, data_folded=True, mask = [True] + [False]*len(sfs) + [True]*addZeroLength)
    else:
        if sampleSize != len(sfs) + 1: return("sampleSize not equal 1+length(sfs)!")
        sfsout = [0.] + sfs + [0.]
        if maskSingletons==True:
            sfsout = dadi.Spectrum(sfsout, data_folded=False, mask = [True]*2 + [False]*(len(sfs)-1))
        else:
            sfsout = dadi.Spectrum(sfsout, data_folded=False, mask = [True] + [False]*len(sfs))
    return sfsout

# One epoch

def one_epoch((nu, T), ns, pts):    
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)    	
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

one_epoch_ex = dadi.Numerics.make_extrap_log_func(one_epoch)
  
def one_epoch_selection((gamma, h), ns, pts):  
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)  	
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))  
    fs[fs < 0] = -fs[fs < 0]
    return fs
    
one_epoch_selection_ex = dadi.Numerics.make_extrap_func(one_epoch_selection, fail_mag=10, extrap_log=True)

# Two epoch

def two_epoch((nu, T), ns, pts):    
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)    
    phi = dadi.Integration.one_pop(phi, xx, T, nu)	
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

two_epoch_ex = dadi.Numerics.make_extrap_log_func(two_epoch)
  
def two_epoch_selection((nu, T, gamma, h), ns, pts):  
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)  
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma, h=h)	
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))  
    fs[fs < 0] = -fs[fs < 0]
    return fs
    
two_epoch_selection_ex = dadi.Numerics.make_extrap_func(two_epoch_selection, fail_mag=10, extrap_log=True)

# Han Danes

def han_danes((N1, N2, T2, NC, TC), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, 0.01, N1)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2)
    dadi.Integration.timescale_factor = 0.000001
    nu_func = lambda t: N2*numpy.exp(numpy.log(NC/N2)*t/TC)
    phi = dadi.Integration.one_pop(phi, xx, TC, nu_func)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
 
han_danes_ex = dadi.Numerics.make_extrap_log_func(han_danes)
    
def han_danes_selection((N1, N2, T2, NC, TC, gamma, h), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, 0.01, N1, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2, gamma=gamma, h=h)
    dadi.Integration.timescale_factor = 0.000001
    nu_func = lambda t: N2*numpy.exp(numpy.log(NC/N2)*t/TC)
    phi = dadi.Integration.one_pop(phi, xx, TC, nu_func, gamma=gamma, h=h)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))   
    fs[fs < 0] = -fs[fs < 0]
    return fs
    
# Three epoch
    
def three_epoch((N1, T1, N2, T2), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, N1)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
 
three_epoch_ex = dadi.Numerics.make_extrap_log_func(three_epoch)
    
def three_epoch_selection((N1, T1, N2, T2, gamma, h), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T1, N1, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2, gamma=gamma, h=h)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    fs[fs < 0] = -fs[fs < 0]
    return fs

# Four epoch
    
def four_epoch((N1, T1, N2, T2, N3, T3), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, N1)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2)
    phi = dadi.Integration.one_pop(phi, xx, T3, N3)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
 
four_epoch_ex = dadi.Numerics.make_extrap_log_func(four_epoch)
    
def four_epoch_selection((N1, T1, N2, T2, N3, T3, gamma, h), ns, pts):
    dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T1, N1, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T2, N2, gamma=gamma, h=h)
    phi = dadi.Integration.one_pop(phi, xx, T3, N3, gamma=gamma, h=h)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    fs[fs < 0] = -fs[fs < 0]
    return fs
    
    
def selection_DFE_gridsearch(sfs_NS, theta, params, ns, demo_selection_dist, sel_dist, output_file):
    
    plist = []
    pcomb = []
    for i in params:
        plist.append(numpy.linspace(i[0], i[1], i[2]).tolist())

    pcomb = itertools.product(*plist)
    
    out = []
    for p in pcomb:

        model = demo_selection_dist(p, ns, sel_dist, theta)
        #print model
        llpois=dadi.Inference.ll(model, sfs_NS)
        llmult=dadi.Inference.ll_multinom(model,sfs_NS)
        #print llpois, llmult

        if not math.isnan(float(llpois)):
            out.append(list(p) + [llpois, llmult])
    
    if output_file:
        fileOut = open(output_file, 'w')
        for i in out:
            print >> fileOut, '{0}\t{1}\t{2}\t{3}'.format(i[0], i[1], i[2], i[3])
        fileOut.close()
    else:
        return(out)
    


