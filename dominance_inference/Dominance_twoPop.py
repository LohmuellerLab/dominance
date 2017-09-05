"""
Functions for two-population inference of dominance
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

try:
    import cPickle as pickle
except:
    import pickle

import types

def copy_func(f, defaults=None):
    defaults = (defaults,)
    return types.FunctionType(f.func_code, f.func_globals, f.func_name, defaults or f.func_defaults, f.func_closure)


#: Counts calls to object_func
_counter = 0
#: Returned when object_func is passed out-of-bounds params or gets a NaN ll.
_out_of_bounds_val = -1e8

# Assume first population (data1) is the outcrossing one, and second population (data2) is the selfing one
# Idea, assume h=0.5 and s twice as large in inbreeding vs. outcrossing population
# optimized are the DFE of s and the s-h relationship (h_intercept and h_rate)
# Note that everything is on scale of s, not Nes

def _object_func(params, data1, data2, cache1, cache2, model_func, sel_dist, scal_fac1, scal_fac2, theta1, theta2,
                 lower_bound=None, upper_bound=None, 
                 verbose=0, multinom=False, flush_delay=0,
                 func_args=[], func_kwargs={}, fixed_params1=None, fixed_params2=None, ll_scale=1,
                 output_stream=sys.stdout, store_thetas=False):
    """
    Objective function for optimization.
    """
    global _counter
    _counter += 1
    
    # Scaling factors scales sel_dist differently for species 1 and species 2

    sel_dist1 = copy_func(sel_dist, defaults = scal_fac1)   # scal_fac1 should be 2*Nea of pop 1
    sel_dist2 = copy_func(sel_dist, defaults = scal_fac2)   # scal_fac2 should be 4*Nea of pop 2
    
    # Deal with fixed parameters
    params_up1 = Inference._project_params_up(params, fixed_params1)
    params_up2 = Inference._project_params_up(params, fixed_params2)
    
    # Check our parameter bounds
    if lower_bound is not None:
        for pval,bound in zip(params_up1, lower_bound):
            if bound is not None and pval < bound:
                return -_out_of_bounds_val/ll_scale
    if upper_bound is not None:
        for pval,bound in zip(params_up1, upper_bound):
            if bound is not None and pval > bound:
                return -_out_of_bounds_val/ll_scale

    ns1 = data1.sample_sizes
    ns2 = data2.sample_sizes
    all_args1 = [params_up1, ns1, sel_dist1, theta1, cache1] + list(func_args)
    all_args2 = [params_up2, ns2, sel_dist2, theta2, cache2] + list(func_args)
    # Pass the pts argument via keyword, but don't alter the passed-in 
    # func_kwargs
    #func_kwargs = func_kwargs.copy()
    #func_kwargs['pts'] = pts
    sfs1 = model_func(*all_args1, **func_kwargs)
    sfs2 = model_func(*all_args2, **func_kwargs)
    if multinom:
        result = Inference.ll_multinom(sfs1, data1) + Inference.ll_multinom(sfs2, data2)
    else:
        result = Inference.ll(sfs1, data1) + Inference.ll(sfs2, data2)
    
    # Bad result
    if numpy.isnan(result):
        result = _out_of_bounds_val
    
    if (verbose > 0) and (_counter % verbose == 0):
        param_str = 'array([%s])' % (', '.join(['%- 12g'%v for v in params_up1]))
        output_stream.write('%-8i, %-12g, %s%s' % (_counter, result, param_str,
                                                   os.linesep))
        Misc.delayed_flush(delay=flush_delay)
    
    return -result/ll_scale

def optimize_log(p0, data1, data2, cache1, cache2, model_func, sel_dist, scal_fac1, scal_fac2 , theta1, theta2, lower_bound=None, upper_bound=None,
                 verbose=0, flush_delay=0.5, epsilon=1e-3, 
                 gtol=1e-5, multinom=False, maxiter=None, full_output=False,
                 func_args=[], func_kwargs={}, fixed_params1=None, fixed_params2=None, ll_scale=1,
                 output_file=None, amoeba=False):
    
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout
    
    args = (data1, data2, cache1, cache2, model_func, sel_dist, scal_fac1, scal_fac2, theta1, theta2, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params1, fixed_params2,
            ll_scale, output_stream)
    
    p0 = Inference._project_params_down(p0, fixed_params1)
    
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
    
    xopt = Inference._project_params_up(numpy.exp(xopt), fixed_params1)
    
    if output_file:
        output_stream.close()
    
    if not full_output or amoeba == True:
        return [-fopt, xopt]
    else:
        return xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag


# Not yet adapted:
def optimize(p0, data1, data2, cache1, cache2, model_func, sel_dist, scal_fac1, scal_fac2 , theta1, theta2, lower_bound=None, upper_bound=None,
                 verbose=0, flush_delay=0.5, epsilon=1e-3, 
                 gtol=1e-5, multinom=False, maxiter=None, full_output=False,
                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
                 output_file=None):
    
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout
    
    args = (data1, data2, cache1, cache2, model_func, sel_dist, scal_fac1, scal_fac2, theta1, theta2, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params, 
            ll_scale, output_stream)
    
    p0 = Inference._project_params_down(p0, fixed_params)
    outputs = scipy.optimize.fmin_bfgs(_object_func, 
                                       p0, epsilon=epsilon,
                                       args = args, gtol=gtol, 
                                       full_output=True,
                                       disp=False,
                                       maxiter=maxiter)
    xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag = outputs
    xopt = Inference._project_params_up(xopt, fixed_params)
    
    if output_file:
        output_stream.close()
    
    if not full_output:
        return [-fopt, xopt]
    else:
        return xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag
        

def _object_func_log(log_params, *args, **kwargs):
    """
    Objective function for optimization in log(params).
    """
    return _object_func(numpy.exp(log_params), *args, **kwargs)

##end of dadi.Inference code



 

