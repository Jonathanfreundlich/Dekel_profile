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
    fun = np.log10(prf.rho_function(r, (c,a,0.5,3,Rvir,Mvir),model='dekel'))
        
    if (a+3.5*np.sqrt(c*10**(-2.))<0.) or (c<0.):
        penalization=10000.
    else:
        penalization=0.
        
    return w*(fun - logrho-penalization)

def errfunc_gnfw(p, r, logrho,Rvir,Mvir,w=1):
    c,a=p 
    fun = np.log10(prf.rho_function(r, (c,a,Rvir,Mvir),model='enfw'))
        
    return w*(fun - logrho)
    
def get_sigma_model(y,ymodel):
    residuals=(y-ymodel)
    N=size(y)
    return np.sqrt(sum(residuals**2)/N)
    
# DEKEL RHO FIT
pfit, perr = optimize.leastsq(errfunc_dekel, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
c_dekel,a_dekel=pfit 
params_dekel_rho=(c_dekel,a_dekel,2,3.,Rvir,Mvir) 
print 'dekel_fit  :', params_dekel_rho

# EINASTO FIT
pfit, perr = optimize.leastsq(errfunc_einasto, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
c_einasto,n_einasto=pfit
params_einasto=(c_einasto,n_einasto,Rvir,Mvir)
print 'einasto_fit :', params_einasto

# ENFW FIT
pfit, perr = optimize.leastsq(errfunc_gnfw, (10.,1.), args=(r,np.log10(rho),Rvir,Mvir,1.), full_output=0)
c_gnfw,a_gnfw=pfit
params_enfw=(c_gnfw,a_gnfw,rvir,mvir) 
print 'gnfw_fit :', params_gnfw

