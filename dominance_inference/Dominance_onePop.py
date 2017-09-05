"""
Estimation of h using dadi.
Assume:h
1) Functional relation between h and s of the form 1/(1/h_intercept + h_rate*s)
2) 2D distribution of s and h?
"""

import os, sys
import numpy
from numpy import logical_and, logical_not
from scipy.special import gammaln
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize
import functools
from dadi import Numerics, Inference, Misc
from dadi.Spectrum_mod import Spectrum
from scipy.interpolate import CubicSpline
import mpmath
import math
from joblib import Parallel, delayed
import multiprocessing
import shelve
import gc

try:
    import cPickle as pickle
except:
	import pickle

#: Stores thetas
_theta_store = {}
#: Counts calls to object_func
_counter = 0
#: Returned when object_func is passed out-of-bounds params or gets a NaN ll.
_out_of_bounds_val = -1e8


# Object function
# Returns log likelihood for a given set of parameters of the DFE and h_intercept and h_rate

def _object_func(params, data, model_func, sel_dist, theta,
                 lower_bound=None, upper_bound=None, 
                 verbose=0, multinom=False, flush_delay=0,
                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
                 output_stream=sys.stdout, store_thetas=False):
    """
    Objective function for optimization.
    """
    global _counter
    _counter += 1
    
    # Deal with fixed parameters
    params_up = Inference._project_params_up(params, fixed_params)
    
    # Check our parameter bounds
    if lower_bound is not None:
        for pval,bound in zip(params_up, lower_bound):
            if bound is not None and pval < bound:
                return -_out_of_bounds_val/ll_scale
    if upper_bound is not None:
        for pval,bound in zip(params_up, upper_bound):
            if bound is not None and pval > bound:
                return -_out_of_bounds_val/ll_scale
    
    ns = data.sample_sizes 
    all_args = [params_up, ns, sel_dist, theta] + list(func_args)

    sfs = model_func(*all_args, **func_kwargs)
    if multinom:
        result = Inference.ll_multinom(sfs, data)
    else:
        result = Inference.ll(sfs, data)
    
    if store_thetas:
        global _theta_store
        _theta_store[tuple(params)] = optimal_sfs_scaling(sfs, data)

    # Bad result
    if numpy.isnan(result):
        result = _out_of_bounds_val
    
    if (verbose > 0) and (_counter % verbose == 0):
        param_str = 'array([%s])' % (', '.join(['%- 12g'%v for v in params_up]))
        output_stream.write('%-8i, %-12g, %s%s' % (_counter, result, param_str,
                                                   os.linesep))
        Misc.delayed_flush(delay=flush_delay)
    
    return -result/ll_scale

# Log version

def _object_func_log(log_params, *args, **kwargs):
    """
        Objective function for optimization in log(params).
        """
    return _object_func(numpy.exp(log_params), *args, **kwargs)

# Optimization function:
# Takes object function to optimize parameters

def optimize_log(p0, data, model_func, sel_dist, theta, lower_bound=None, upper_bound=None,
                 verbose=0, flush_delay=0.5, epsilon=1e-3, 
                 gtol=1e-5, multinom=False, maxiter=None, full_output=False,
                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
                 output_file=None, amoeba = False):
    
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout
    
    args = (data, model_func, sel_dist, theta, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params, 
            ll_scale, output_stream)
    
    p0 = Inference._project_params_down(p0, fixed_params)
    
    if amoeba == False:
        
        outputs = scipy.optimize.fmin_bfgs(_object_func_log, 
                                           numpy.log(p0), epsilon=epsilon,
                                           args = args, gtol=gtol, 
                                           full_output=True,
                                           disp=False,
                                           maxiter=maxiter)
        xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag = outputs
        
    else:

        outputs = scipy.optimize.fmin(_object_func_log, 
                                           numpy.log(p0),
                                           args = args, 
                                           full_output=True,
                                           disp=False,
                                           maxiter=maxiter)
        xopt, fopt, iter, funcalls, warnflag = outputs
        
    xopt = Inference._project_params_up(numpy.exp(xopt), fixed_params)
    
    if output_file:
        output_stream.close()
    
    if not full_output or amoeba == True:
        return [-fopt, xopt]
    else:
        return xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag

# Function for finding nearest values in array:
        
def find_nearest_idx(array,value):
    idx = numpy.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx      
        
# Returns SFS given a DFE parameters, sample size, a DFE distribution, and theta.
# Basically, this integrates across Nes*h, params are first DFE then h_intercept and h_rate

def demo_selection_dist(params, ns, sel_dist, theta, cache):
    """
    sel_dist should be a function that is evaluated PDF(X) = func(x, param1, param2 ..., paramn)
    
    theta is required for now, just going to use Poisson ll
    """
    #load saved objects
    
    #spectra_obj = pickle.load(open('{0}spectra.obj'.format(note),'rb'))
        
    spectra_obj = cache
    
    # Note that first and last entry of SFS are meaningless! 
    # The last two entries of params are now h_intercept and h_rate
    
    params_DFE = params[:-2]
    params_h = params[-2:]

    hvalues = comp_h_from_s(spectra_obj['gammas'], *params_h)
      
    # Choose the closest SFS that correspond to the respective hvalues:
        
    hlist_idx = [find_nearest_idx(spectra_obj['hlist'], h) for h in hvalues]
    
    spectra_interp = numpy.array([spectra_obj['spectra'][hlist_idx[i], i, :] for i in range(spectra_obj['spectra'].shape[1])])
                 
    # For some reason, some SFS just contain nan when 2Neas*h is large...
    # For now, set them to 0 and deal with it later... they should only contain small values anyway
    # Update: THis is fixed, the problem where negative values in the SFS
    
    # spectra_interp = [numpy.nan_to_num(sfs) for sfs in spectra_interp]  # replaces nan with 0
             
    #compute weights for each SFS
    sel_args = (spectra_obj['gammas'],) + tuple(params_DFE)
    weights = sel_dist(*sel_args)
    
    #compute weight for the effectively neutral portion. not using CDF function because
    #I want this to be able to compute weight for an arbitrary mass functions
    weight_neu, err_neu = scipy.integrate.quad(sel_dist, spectra_obj['gammas'][-1], 0, args=tuple(params_DFE))
    weight_lethal, err_lethal = scipy.integrate.quad(sel_dist, -numpy.inf, spectra_obj['gammas'][0], args=tuple(params_DFE))
    
    #function's adaptable for demographic models from 1-3 populations
    pops = len(spectra_obj['neu_spec'].shape)
    if pops == 1:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    elif pops == 2:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    elif pops == 3:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis,numpy.newaxis,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    else:
        raise IndexError("Must have one to three populations")
    
    integrated_fs = Spectrum(integrated, extrap_x=spectra_obj['extrap_x'])

    # Changed this:
    # Lethal mutations now don't fall out. All lethal mutations contribute the most deleterious SFS.
    # This assumes that the range of gamma goes from small to lethal!
    
    return integrated_fs * theta


# Here, for now the above function is just copied...
def demo_selection_distINV(params, ns, sel_dist, theta, cache):
    """
    sel_dist should be a function that is evaluated PDF(X) = func(x, param1, param2 ..., paramn)
    
    theta is required for now, just going to use Poisson ll
    """
    #load saved objects
    
    #spectra_obj = pickle.load(open('{0}spectra.obj'.format(note),'rb'))
        
    spectra_obj = cache
    
    # Note that first and last entry of SFS are meaningless! 
    # The last two entries of params are now h_intercept and h_rate
    
    params_DFE = params[:-2]
    params_h = params[-2:]

    hvalues = comp_h_from_s_INV(spectra_obj['gammas'], *params_h)
      
    # Choose the closest SFS that correspond to the respective hvalues:
        
    hlist_idx = [find_nearest_idx(spectra_obj['hlist'], h) for h in hvalues]
    
    spectra_interp = numpy.array([spectra_obj['spectra'][hlist_idx[i], i, :] for i in range(spectra_obj['spectra'].shape[1])])
                 
    # For some reason, some SFS just contain nan when 2Neas*h is large...
    # For now, set them to 0 and deal with it later... they should only contain small values anyway
    # Update: THis is fixed, the problem where negative values in the SFS
    
    # spectra_interp = [numpy.nan_to_num(sfs) for sfs in spectra_interp]  # replaces nan with 0
             
    #compute weights for each SFS
    sel_args = (spectra_obj['gammas'],) + tuple(params_DFE)
    weights = sel_dist(*sel_args)
    
    #compute weight for the effectively neutral portion. not using CDF function because
    #I want this to be able to compute weight for an arbitrary mass functions
    weight_neu, err_neu = scipy.integrate.quad(sel_dist, spectra_obj['gammas'][-1], 0, args=tuple(params_DFE))
    weight_lethal, err_lethal = scipy.integrate.quad(sel_dist, -numpy.inf, spectra_obj['gammas'][0], args=tuple(params_DFE))
    
    #function's adaptable for demographic models from 1-3 populations
    pops = len(spectra_obj['neu_spec'].shape)
    if pops == 1:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    elif pops == 2:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    elif pops == 3:
        integrated = spectra_obj['neu_spec']*weight_neu + Numerics.trapz(weights[:,numpy.newaxis,numpy.newaxis,numpy.newaxis]*spectra_interp, spectra_obj['gammas'], axis=0) + spectra_interp[0]*weight_lethal
    else:
        raise IndexError("Must have one to three populations")
    
    integrated_fs = Spectrum(integrated, extrap_x=spectra_obj['extrap_x'])

    # Changed this:
    # Lethal mutations now don't fall out. All lethal mutations contribute the most deleterious SFS.
    # This assumes that the range of gamma goes from small to lethal!
    
    return integrated_fs * theta
    
    
# Generate SFS for a grid of values of Neas and h
# For now, generate SFS for 21 h values from 0 to 1

def interp_sfs_freq(gammas, freq, gammas_ex, logx=True):
    freq = numpy.array(freq)
    # Some very small negative values appear and make problems. Make them positive for now:
    # freq[freq < -1e-50] = -freq[freq < -1e-50]
    freq[freq < 0] = -freq[freq < 0]
    # Cubic interpolation in log-log-space:
    if logx == True:
        fint = CubicSpline(-numpy.log(-gammas), numpy.log(freq))
        return [gammas_ex, numpy.exp(fint(-numpy.log(-gammas_ex)))]
    # Cubic interpolation in linear-log-space:
    else:
        fint = CubicSpline(gammas, numpy.log(freq))
        return [gammas_ex, numpy.exp(fint(gammas_ex))]

def interp_full_sfs(gammas, sfss, gammas_ex, logx=True):  # I don't use this function, but might be useful (it was too slow)
    sfss = numpy.array(sfss)
    # Cubic interpolation:
    if logx == True:
        outsfss = numpy.array([interp_sfs_freq(gammas, freq, gammas_ex, logx=True)[1] for freq in sfss.transpose()])
        return [gammas_ex, outsfss]
    else:
        outsfss = numpy.array([interp_sfs_freq(gammas, freq, gammas_ex, logx=False)[1] for freq in sfss.transpose()])
        return [gammas_ex, outsfss]

def interp_spectra_gamma(spectra, gammas, New_Npts):
    gammas_ex = -numpy.logspace(numpy.log10(-gammas[-1]), numpy.log10(-gammas[0]), New_Npts)[::-1]
    spectra_new = numpy.ndarray(shape = [spectra.shape[0], New_Npts, spectra.shape[2]], dtype=numpy.float32)
    for h_index in range(spectra.shape[0]):
        for daf in range(1, spectra.shape[2]-1):
            freq = [spectra[h_index][i][daf] for i in range(spectra.shape[1])]
            interp_freq = interp_sfs_freq(gammas, freq, gammas_ex, logx=True)
            spectra_new[h_index,:,daf] = interp_freq[1]
    return([gammas_ex, spectra_new])
    
def interp_spectra_h(spectra, hlist, New_Npts):
    hlist_ex = numpy.linspace(hlist[0], hlist[-1], New_Npts)
    spectra_new_new = numpy.ndarray(shape = [New_Npts, spectra.shape[1], spectra.shape[2]], dtype=numpy.float32)
    for Ns_index in range(spectra.shape[1]):
        for daf in range(1, spectra.shape[2]-1):
            freq = [spectra[i][Ns_index][daf] for i in range(spectra.shape[0])]
            interp_freq = interp_sfs_freq(hlist, freq, hlist_ex, logx=False)
            spectra_new_new[:,Ns_index,daf] = interp_freq[1]
    return([hlist_ex, spectra_new_new])
    
    
def processInput(params, gammas, h, ns, pts, demo_sel_func, neu_singFreq):
    spectra_temp = [None] * len(gammas)
    flag_strong_del = False
    for ii,gamma in enumerate(gammas[::-1]):
        if flag_strong_del == False:
            sfs = demo_sel_func(params+(gamma, h), ns, pts)
            spectra_temp[ii] = sfs
            if sfs[1] < neu_singFreq*1e-4:  # When singletons are less than 1e-4*neutral frequency then assume SFS don't contribute
                flag_strong_del = True
                strong_del_sfs = sfs*0 + 1e-30
        else:
            spectra_temp[ii] = strong_del_sfs
        
        # Print out some info:
        print "Finished: "+str(h)+" "+str(gamma)+" Ignore singletons? "+str(flag_strong_del)  # Remove this         
        if any(numpy.isnan(spectra_temp[ii].flatten())):
            print "Computed sfs contains nan, I will set those entries to 1e-30"
        if any(numpy.isinf(spectra_temp[ii].flatten())):
            print "Computed sfs contains inf, I will set those entries to 1e-30"
            
        sys.stdout.flush()   # Remove this 
    return spectra_temp[::-1]

def generatecache(params, ns, demo_sel_func, note="", pts=1000, Npts_large = 1000, Npts_h_large = 1001, NeAnc=7000, num_cores = 1):
    """
    params need to be demographic params without gamma
    spectra are pickled and saved into working directory
    better pickling management to come
    """
    
    print("Note that the upper integration bound is the eff. population size (homozygous lethal mutations)!")
    
    int_bounds = tuple([1e-4, NeAnc])
    
    largestGamma = int_bounds[1]; smallestGamma = int_bounds[0]
    upperExponentBound = numpy.log10(largestGamma)
    lowerExponentBound = numpy.log10(smallestGamma)
    Npts_small = numpy.ceil(5*(upperExponentBound-lowerExponentBound))   # 5 points per decimal exponent works well. Hardcoded
    
    gammas = -numpy.logspace(numpy.log10(int_bounds[1]), numpy.log10(int_bounds[0]), Npts_small)
    hlist = numpy.linspace(0.0, 1.0, num=21)  # 20 points of h between 0 and 1 works well. Hardcoded

    #hlist[numpy.where(hlist == 0.5)] = 0.5001 # Make sure h is not exactly 0.5. Seems that additive model does not work in vanilla dadi for some reason when selection strong

    #generate a neutral fs
    neu_spec = numpy.array(demo_sel_func(params+(0,0.5), ns, pts))
   
    neu_singFreq = neu_spec[1]
    #demo_ex = dadi.Numerics.make_extrap_log_func(demo_func)
    	
    # make spectra, a two-dimensional list of SFSs with different h and 2Neas
    # outer dimension: h goes from recessive (0) to dominant (1)
    # inner dimension: 2Neas goes from deleterious to neutral
    
    spectra = []

    # This uses multiprocessing and also ignores strongly deleterious SFSs:
    
    # num_cores = multiprocessing.cpu_count()
    spectra = Parallel(n_jobs=num_cores)(delayed(processInput)(params, gammas, h, ns, pts, demo_sel_func, neu_singFreq) for h in hlist)
    
    
    extrap_x = spectra[0][0].extrap_x

    spectra = numpy.array(spectra, dtype=numpy.float32)
    
    spectra_obj = dict(extrap_x = extrap_x, spectra = spectra, neu_spec = neu_spec, gammas = gammas, hlist = hlist)
    print "Writing out not-yet-interpolated file..."
    out = open('{0}spectra_notInterpolated.obj'.format(note),'wb')
    pickle.dump(spectra_obj, out)
    out.close()
    
    neu_spec, spectra, gammas, hlist, extrap_x = pickle.load(open('{0}spectra_notInterpolated.obj'.format(note),'rb')).values()
    
    spectra[spectra <= 1e-30] = 1e-30   # Entries can be zero and make problems with interploation in log space. Set them to 1e-30.
    spectra[numpy.isnan(spectra)] = 1e-30
    spectra[numpy.isinf(spectra)] = 1e-30
            
           
    # Interpolate spectra to get Npts_large number of gamma SFS's
    
    spectra_new = interp_spectra_gamma(spectra, gammas, Npts_large)
    
    # Interpolate spectra to get 1001 number of h SFS's (hardcoded for now)
    
    if any(numpy.isnan(spectra_new[1]).flatten()):
        print "Computed spectra_new contain nan, I set those entries to 1e-30"
        spectra_new[1][numpy.isnan(spectra_new[1])] = 1e-30   # There were some problems with strong selection leading to nan... set to 1e-30.
    
    if any(numpy.isinf(spectra_new[1]).flatten()):
        print "Computed spectra_new contain inf, I set those entries to 1e-30"
        spectra_new[1][numpy.isinf(spectra_new[1])] = 1e-30
    
    spectra_new_new = interp_spectra_h(spectra_new[1], hlist, Npts_h_large)    # Old code: working, but takes up too much memory for large SFS. Solution: Dont pickle large array, and use dtype32 instead of dtype64!!

    # Here, set all the strange sfs entries that are from strong selection (NS<-20000) AND dominance (h>0.3) to 1e-30
    spectra_new_new[1][numpy.where(spectra_new_new[0] > 0.3)[0][0]:, :numpy.where(spectra_new[0] < -20000)[0][-1], :] = 1e-30
 
    #save info to avoid having to recompute each time a different distribution is evaluated
    
    spectra_obj = dict(extrap_x = extrap_x, neu_spec = neu_spec, arrayFN = '{0}spectra_bigArray.npy'.format(note), gammas = spectra_new[0], hlist = spectra_new_new[0])

    print "Writing out final files..."
    out = open('{0}spectra_info.obj'.format(note),'wb')
    pickle.dump(spectra_obj, out)
    out.close()
    
    numpy.save('{0}spectra_bigArray.npy'.format(note), spectra_new_new[1])


def loadSpectraObj(note):
    spectra_obj = pickle.load(open('{0}spectra_info.obj'.format(note),'rb'))
    spectra_obj['spectra'] = numpy.load(spectra_obj['arrayFN'])
    return spectra_obj


#############################
# h-s relationship function #
#############################

# Assume a specific h-s relationship
# 0 < h_intercept <= 1
# 0 <= h_rate <= 100,000

# Works also if s is Neas instead

def comp_h_from_s(s, h_intercept, h_rate):
    if h_intercept==0:
        return([0 for i in s])
    hout = 1/(1/h_intercept + h_rate*abs(s))
    return(hout)

    
def comp_h_from_s_INV(s, h_intercept, h_rate):   # Function for increasing h with more deleterious s, converging to h=1
    if h_intercept==1:
        return([1 for i in s])
    hout = 1-1/(1/(1-h_intercept) + h_rate*abs(s))
    return(hout)
    
#############################
# Default DFE distributions #
#############################


def neutMass_func(mgamma, cutoff=(-1e-4)):
    if mgamma > cutoff and mgamma <= 0:   
        return -1/cutoff
    else:
        return 0

def neutMass(mgamma, scal_fac=1):
    if isinstance(mgamma, int) or isinstance(mgamma, float):
        return(neutMass_func(mgamma, -1e-4*scal_fac))
    return(numpy.array([neutMass_func(gamma, -1e-4*scal_fac) for gamma in mgamma]))
            

def gamma_dist(mgamma, alpha, beta, scal_fac=1):
    """
    x, shape, scale
    """
    return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta*scal_fac)
    
def shifted_gamma_dist(mgamma, alpha, beta, shift, scal_fac=1):   
    """
    x, shape, scale
    """
    return((scipy.stats.distributions.gamma.pdf((-mgamma + shift)*scal_fac, alpha, scale=beta))/scal_fac)
    
def gamma_neutMass_dist(mgamma, alpha, beta, pneut, scal_fac=1):
    """
    x, shape, scale
    """ 
    return((1-pneut) * scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta*scal_fac) + pneut * neutMass(mgamma, scal_fac=scal_fac))
    
def beta_dist(mgamma, alpha, beta):
    """
    x, alpha, beta
    """
    return scipy.stats.distributions.beta.pdf(-mgamma, alpha, beta)

def exponential_dist(mgamma, scale):
    return scipy.stats.distributions.expon.pdf(-mgamma, scale=scale)

def lognormal_dist(mgamma, mu, sigma, scal_fac=1):
    return scipy.stats.distributions.lognorm.pdf(-mgamma, sigma, scale=numpy.exp(mu + numpy.log(scal_fac)))

def normal_dist(mgamma, mu, sigma):
    return scipy.stats.distributions.norm.pdf(-mgamma, loc=mu, scale=sigma)
 
 
from mpmath import mp
mp.dps = 50 
 
def gamma_rice_dist(mgamma, alpha, beta, NeInv, Ne_scal, scal_fac=1):  # Note that Ne_scal should be Ne_anc
    if isinstance(mgamma, int) or isinstance(mgamma, float):
        s = mgamma/2./Ne_scal*scal_fac
        prob = gamma_dist(-numpy.abs(s), alpha, beta, scal_fac=1)/(1.+mp.exp(2/NeInv*s))
        return prob/2./Ne_scal*scal_fac
    s = [gamma/2./Ne_scal*scal_fac for gamma in mgamma]
    prob = [gamma_dist(-numpy.abs(si), alpha, beta, scal_fac=1)/(1.+mp.exp(2/NeInv*si)) for si in s]
    return numpy.array([p/2./Ne_scal*scal_fac for p in prob])
    
    
    
 # The Lourenco et al., 2011 function using mpmath

# In equilibrium
  

def lourenco_eq_dist_mp(mgamma, m, sigma, Ne, Ne_scal, scal_fac=1):   # not sure if I do it right with the scale factor...
     s = mgamma/2./Ne_scal*scal_fac
     #s = mgamma
     prob = mp.power(2, (1-m)/2.)*mp.power(Ne, 0.5)*mp.power(mp.fabs(s), (m-1)/2.)*(1+1/(Ne*mp.power(sigma, 2.)))*mp.exp(-Ne*s) / (mp.power(mp.pi,0.5)*mp.power(sigma, m)*mp.gamma(m/2.)) * mp.besselk((m-1)/2., Ne*mp.fabs(s)*mp.power(1+1/(Ne*mp.power(sigma, 2.)),0.5))
     return float(prob/2/Ne_scal*scal_fac)
     #return float(prob)

def lourenco_eq_dist(mgamma, m, sigma, Ne, Ne_scal, scal_fac=1):
    if isinstance(mgamma, int) or isinstance(mgamma, float):
        return lourenco_eq_dist_mp(mgamma, m, sigma, Ne, Ne_scal, scal_fac)
    out = [lourenco_eq_dist_mp(gamma, m, sigma, Ne, Ne_scal, scal_fac) for gamma in mgamma]
    return numpy.array(out)
    
# Not in equilibrium
    
def lourenco_dist_mp(mgamma, z, n, m, sigma, Ne, Ne_scal, scal_fac=1):   # not sure if I do it right with the scale factor...
     s = mgamma/2./Ne_scal*scal_fac
     if s==0:
         return numpy.nan
     #s = mgamma
     # prob = mp.power(2, (1-m)/2.)*mp.power(Ne, 0.5)*mp.power(mp.fabs(s), (m-1)/2.)*(1+1/(Ne*mp.power(sigma, 2.)))*mp.exp(-Ne*s) / (mp.power(mp.pi,0.5)*mp.power(sigma, m)*mp.gamma(m/2.)) * mp.besselk((m-1)/2., Ne*mp.fabs(s)*mp.power(1+1/(Ne*mp.power(sigma, 2.)),0.5))
     prob = mp.power(2, -(m+1)/2)*mp.exp(-n*s/4/mp.power(z, 2.))*mp.sqrt(n)*mp.power((1+4.*mp.power(z, 2.)/n/mp.power(sigma, 2.))/mp.power(s, 2.), (1-m)/4)*mp.power(sigma, -m) / (mp.sqrt(mp.pi)*z*mp.gamma(m/2.)) * mp.besselk((m-1)/2., 1/(4.*mp.sqrt(mp.power(z, 2.)/(n*mp.power(s,2.)*(n/mp.power(z,2.)+4./mp.power(sigma,2.))))))

     return float(prob/2/Ne_scal*scal_fac)
     #return float(prob)

def lourenco_dist(mgamma, z, n, m, sigma, Ne, Ne_scal, scal_fac=1):
    if isinstance(mgamma, int) or isinstance(mgamma, float):
        return lourenco_dist_mp(mgamma, z, n, m, sigma, Ne, Ne_scal, scal_fac)
    out = [lourenco_dist_mp(gamma, z, n, m, sigma, Ne, Ne_scal, scal_fac) for gamma in mgamma]
    return numpy.array(out)
    

