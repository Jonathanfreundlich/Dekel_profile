###################################################################
#                    FITTING THE DEKEL+ PROFILE                   #
###################################################################

#      This python program enables to fit density profiles        #
#      with the Dekel+ profile and other profiles                 # 
#              Written by J. Freundlich, 22/04/2020               #
#                                                                 #
#      Note: we impose a negative slope at the resolution         #
#      for the Dekel+ profile                                     #
#      - r: radius (array)                                        #
#      - rho: density profile (array)                             #

###################################################################

import numpy as np
from scipy import optimize
import Dekel_profile as prf

def errfunc_einasto(p, r, logrho,Rvir,Mvir,w=1):
    c,n=p 
    fun = np.log10(prf.rho_function(r, (c,n,Rvir,Mvir),model='einasto'))      
    return w*(fun - logrho)
    
def errfunc_dekel(p, r, logrho,Rvir,Mvir,w=1):
    c,a=p 
    if (a > 2.999) or (c<0.) or (a+3.5*np.sqrt(c*10**(-2.))<0.):
        return -np.inf * np.ones( len(r) )
    
    fun = np.log10(prf.rho_function(r, (c,a,2,3,Rvir,Mvir),model='dekel'))
        
    return w*(fun - logrho)

def errfunc_gnfw(p, r, logrho,Rvir,Mvir,w=1):
    c,a=p 
    fun = np.log10(prf.rho_function(r, (c,a,Rvir,Mvir),model='gnfw'))
        
    return w*(fun - logrho)

def errfunc_nfw(p, r, logrho,Rvir,Mvir,w=1):
    c=p[0]
    fun = np.log10(prf.rho_function(r, (c,Rvir,Mvir),model='nfw'))
        
    return w*(fun - logrho)
    
def get_sigma_model(y,ymodel):
    residuals=(y-ymodel)
    N=size(y)
    return np.sqrt(sum(residuals**2)/N)
    

# DEKEL rho fit
def best_fit_Dekel( r, rho, Rvir,Mvir ):
    pfit, perr = optimize.leastsq( errfunc_dekel, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
    c_dekel,a_dekel=pfit 
    params_dekel_rho=(c_dekel,a_dekel,2,3.,Rvir,Mvir) 
    return params_dekel_rho

# NFW rho fit
def best_fit_NFW( r, rho, Rvir,Mvir ):
    pfit, perr = optimize.leastsq( errfunc_nfw, (10.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
    c_nfw =  pfit[0]
    params_nfw =  (c_nfw, Rvir,Mvir)
    return params_nfw


# gNFW rho fit
def best_fit_gNFW( r, rho, Rvir,Mvir ):
    pfit, perr = optimize.leastsq( errfunc_gnfw, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
    c_gnfw,a_gnfw=pfit
    params_gnfw=(c_gnfw,a_gnfw,Rvir,Mvir) 
    return params_gnfw


# Einasto rho fit
def best_fit_Einasto( r, rho, Rvir,Mvir ):
    pfit, perr = optimize.leastsq( errfunc_einasto, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
    c_einasto,n_einasto=pfit
    params_einasto=(c_einasto,n_einasto,Rvir,Mvir)
    return params_einasto


