###################################################################
#           ANALYTIC EXPRESSIONS FOR THE DEKEL+ PROFILE           #
###################################################################

#      This python program enables to compute the analytical      #
#      expressions of the Dekel+ dark matter density profile      #
#              Written by J. Freundlich, 22/04/2020               #
#                                                                 #
#      Based on:                                                  #
#      - Zhao (1996), MNRAS, 278, 488                             #
#      - Dekel et al. (2017), MNRAS, 468, 1005                    #
#      - Freundlich et al. (2020), arXiv:2004.08395               #
#                                                                 #
#      params=(c,a,b,g,Mvir,Rvir) with b=2, g=3                   #
#      quantities in units of kpc, Msun, Gyr,...                  #

###################################################################

import numpy as np
from math import factorial
from scipy.special import betainc, gamma
from scipy.integrate import quad

G=4.499753324353496e-06 # [kpc^3 Msun^-1 Gyr^-2]

# DENSITY [Msun kpc^-3]
def rho_function(r, params, model='dekel'):
    
    # DEKEL+ PROFILE
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return (q_bar(c, a, b, g, Mvir, Rvir)*float(3.-a)/3.*(1.+float(3.-g)/(3.-a)*pow(x,1./b))/
                np.array(pow(x,a)*pow((1.+pow(x,1./b)),b*(g-a)+1.)))
    
    # NFW PROFILE
    if model == 'nfw':
        (c, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        qb = Mvir*pow(c,3)/(4.*np.pi*pow(Rvir,3)*(np.log(1.+c)-c/(1.+c)))
        return qb/(x*pow(1.+x,2)) 
    
    # EINASTO PROFILE
    if model == 'einasto':
        (c, n, Rvir, Mvir) = params
        rs = float(Rvir)/c
        h = rs/float(pow(2*n,n))    
        raw_rho = np.exp(-pow((r/h),pow(n,-1)))
        q = Mvir / quad( dM_einasto,1.e-8*Rvir,Rvir,args=(1.,h,n), epsabs=0., epsrel=1.e-5)[0]
        return q*raw_rho
    
    # DI CINTIO PROFILE
    if model == 'dicintio':
        (rho_s,r_s,a,b,g) = params
        x=r/r_s
        return rho_s/(pow(x,g)*pow(1+pow(x,a),(b-g)/a))
    
    # GENERALIZED NFW PROFILE
    if model == 'gnfw':
        (c,a,Rvir,Mvir) = params
        rs=float(Rvir)/c
        x=r/rs
        raw_rho=1./(pow(x,a)*pow(1.+x,3.-a))
        rhoc = Mvir / quad( dM_gNFW,1.e-8*Rvir,Rvir,args=(1.,rs,a), epsabs=0., epsrel=1.e-5)[0]
        return rhoc/(pow(x,a)*pow(1.+x,3.-a))
        
# MEAN DENSITY  [Msun kpc^-3]
def brho_function(r, params, model='dekel'):
    
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return q_bar(c, a, b, g, Mvir, Rvir)/np.array(pow(x,a)*pow((1+pow(x,1./b)),b*(g-a)))
    
    if model == 'nfw':
        return M_function(r, params, model)/(4.*np.pi/3.*pow(r,3))
    
    if model == 'dicintio':
        return M_function(r, params, model)/(4.*np.pi/3.*pow(r,3))
        
# ENCLOSED MASS [Msun]
def M_function(r, params, model='dekel'):
    
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return Mvir*mu(c, a, b, g)*pow(x,3.-a)/np.array(pow(1.+pow(x,1./b),b*(g-a)))
    
    if model == 'nfw':
        (c, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return float(Mvir)*(np.log(1+pow(x,1.))-pow(x,1.)/(1+pow(x,1.)))/(np.log(1.+c)-c/(1.+c))
    
    if model == 'dicintio':
        (rho_s,r_s,a,b,g) = params
        M_array=np.nan*np.ones_like(r)
        for i in range(np.size(r)):
            M_array[i]=quad(dM_dicintio,0,r[i],args=(rho_s,r_s,a,b,g))[0]
        return M_array
   
    if model == 'einasto':
        (c, n, Rvir, Mvir) = params
        rs = float(Rvir)/c
        h = rs/float( pow(2*n,n) )
        q = Mvir / quad( dM_einasto,1.e-8*Rvir,Rvir,args=(1.,h,n), epsabs=0., epsrel=1.e-5)[0]
        M_array=np.nan*np.ones_like(r)
        for i in range(np.size(r)):
            M_array[i]=q * quad( dM_einasto,1.e-8*Rvir,r[i],args=(1.,h,n), epsabs=0., epsrel=1.e-5)[0]
        return M_array
    
    if model == 'gnfw':
        (c,a,Rvir,Mvir) = params
        rs=float(Rvir)/c
        x=r/rs
        rhoc = Mvir / quad( dM_gNFW,1.e-8*Rvir,Rvir,args=(1.,rs,a), epsabs=0., epsrel=1.e-5)[0]
        M_array=np.nan*np.ones_like(r)
        for i in range(np.size(r)):
            M_array[i]= rhoc * quad( dM_gNFW,1.e-8*Rvir,r[i],args=(1.,rs,a), epsabs=0., epsrel=1.e-5)[0]
        return M_array         

# ORBITAL VELOCITY [kpc Gyr^-1]
def V_function(r, params, model='dekel'):
    if model == 'dekel':     
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        Vvir2 = G*Mvir/Rvir
        return np.sqrt(Vvir2*mu(c, a, b, g)*c*pow(x,2.-a)/np.array(pow(1.+pow(x,1./b),b*(g-a))))

# LOGARITHMIC DENSITY SLOPE []
def s_function(r, params, model='dekel'):

    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        return -(3.-g)/float(3.-a)/b*pow(x,1./b)/(1.+(3.-g)/float(3.-a)/b*pow(x,1./b)) + (a+(g+1./b)*pow(x,1./b))/np.array(1.+pow(x,1./b))
        
    if model == 'dicintio':
        (rho_s,r_s,a,b,g) = params
        x = r/r_s
        return g+(b-g)*x**a/(1.+pow(x,a))
    
    if model == 'einasto':
        (c, n, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x=r/rs
        return 2.*pow(x,1./n)
    
    if model == 'gnfw':
        (c,a,Rvir,Mvir) = params
        rs=float(Rvir)/c
        x=r/rs
        return (a+3.*x)/(1+x)
        
# POTENTIAL [kpc^2 Gyr^-2]
def U_function(r, params, model='dekel'): 
    
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        Vvir2 = G*Mvir/Rvir
        
        if (b==1 and g==3):
            return c*mu(c,a,b,g)/float(2.-a)*Vvir2*(pow(x/np.array(1.+x),2.-a)-pow(c/float(1.+c),2.-a))-Vvir2
            
        if (b==2 and g==3): 
            chi = pow(x,0.5)/np.array(1.+pow(x,0.5))
            chic = pow(c,0.5)/np.array(1.+pow(c,0.5))
            return -Vvir2-2*c*mu(c,a,b,g)*Vvir2*((pow(chic,2*(2-a))-pow(chi,2*(2-a)))/float(2*(2-a))-(pow(chic,2*(2-a)+1)-pow(chi,2*(2-a)+1))/float(2*(2-a)+1))
        
        else:
            dU = lambda y: G*M_function(y,params,model)/array(pow(y,2))
            return np.array([-Vvir2-quad(dU,rt,Rvir)[0] for rt in r])

# VELOCITY DISPERSION [kpc Gyr^-1] (ISOTROPIC)
def sigmar_function(r,params, model='dekel'):
    
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        X = chi(x)
        Xc= chi(c)
        factor=2.*G*Mvir/Rvir*c*mu(c, a, b, g)*x**a*pow(1.+pow(x,0.5),2.*(3.5-a))
        u = 4*(1.-a)
        
        if u>0: # slightly more rapid
            sigmar2=factor*(Binc(u,9,1)-Binc(u,9,X))
        else: #(Binc(u,9+2*n,1)-Binc(u,9+2*n,X))
            sigmar2=factor*Dbetainc(u,9,X) 
            
        return np.sqrt(sigmar2)

# VELOCITY DISPERSION [kpc Gyr^-1] (ISOTROPIC) AS A SUM 
# (Zhao 1996, Eq. 19, Freundlich+2020a, Eq. B8)
def sigmar_function_sum(r,params, model='dekel'):

    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        X = chi(x)
        Xc= chi(c)
        Vvir2 = G*Mvir/Rvir
        u = 4*(1.-a)
    
        return np.sqrt(2.*Vvir2*c/pow(Xc,2.*(3.-a))*pow(x,3.5)/pow(X,2.*(3.5-a)) \
            * ( (1.-X**u)/u - 8.*(1.-X**(u+1.))/(u+1.) \
            + 28.*(1.-X**(u+2.))/(u+2.) - 56.*(1.-X**(u+3.))/(u+3.) \
            + 70.*(1.-X**(u+4.))/(u+4.) - 56.*(1.-X**(u+5.))/(u+5.) \
            + 28.*(1.-X**(u+6.))/(u+6.) - 8.*(1.-X**(u+7.))/(u+7.) \
            + (1.-X**(u+8.))/(u+8.)))

# VELOCITY DISPERSION [kpc Gyr^-1] (ISOTROPIC) AS ANOTHER SUM 
# (Freundlich+2020a, Eq. B10)        
def sigmar_function_sum2(r,params, model='dekel'):
      
    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        rs = float(Rvir)/c
        x = r/rs
        X = chi(x)
        Xc= chi(c)
        Vvir2 = G*Mvir/Rvir
        u = 4*(1.-a)
  
        return np.sqrt(2.*Vvir2 *c/pow(Xc,2.*(3.-a)) * x**a*(1.+x**0.5)**(2.*(3.5-a)) \
        * (B9(u,1)-B9(u,X)))
 
# DISTRIBUTION FUNCTION [in Virial units, i.e., Vvir^-3 Rvir^-3]
def distribution_function(E,xE,params, model='dekel'):
    '''
    xE=Psi^-1(Ei,params,model) with Psi=-U the relative potential
    xE can notably be obtained numerically from Psi(xE)=E
    '''

    if model == 'dekel':
        (c, a, b, g, Rvir, Mvir) = params
        if np.size(E)==1:
            fE= 2./(np.sqrt(8.)*np.pi)*quad(distribution_integrand,xE,c, args=(E,params,model))[0]-1/np.sqrt(2)/np.pi*dnudpsi(c,params)/dpsi(c,params)*np.sqrt(E+U_adim(c,params))
        if np.size(E)>1:
            fE=np.nan*np.ones_like(np.array(E))
            for i,Ei in enumerate(E):
                xmin_i=xE[i]
                try:
                    fE[i]=2./(np.sqrt(8.)*np.pi)*quad(distribution_integrand,xmin_i,c, args=(Ei,params,model))[0]-1/np.sqrt(2)/np.pi*dnudpsi(c,params)/dpsi(c,params)*np.sqrt(Ei-psi(c,params))
                except:
                    fE[i]=nan

        return fE
        
###################################################################
# MASS-DEPENDENT PRESCRIPTIONS

# DI CINTIO PROFILE
def get_params_dicintio(rvir,mstar,mvir,fit_rvir='BN'):
    logMM=np.log10(mstar/mvir)
    
    # alpha, beta, gamma from Di Cintio + 14
    a=2.94-np.log10((10**(logMM+2.33))**-1.08+(10**(logMM+2.33))**2.29)
    b=4.23+1.34*logMM+0.26*logMM**2
    g=-0.06+np.log10((10**(logMM+2.56))**-0.68+10**(logMM+2.56))
    
    # DM concentration from Dutton & Maccio 14
    if fit_rvir=='BN':
        c_DM=10**(1.025-0.097*np.log10(mvir/1.e12*0.671)) 
    elif fit_rvir=='R200':
        c_DM=10**(0.905-0.101*np.log10(mvir/1.e12*0.671))
    
    # SPH concentration from Di Cintio + 14
    c_SPH=(1.0+0.00003* np.exp(3.4*(logMM+4.5)))*c_DM
    
    r_2=rvir/c_SPH
    r_s=r_2*((b-2.)/(2.-g))**(1./a)

    rho_s=mvir/quad(dM_dicintio,0,rvir,args=(1.,r_s,a,b,g))[0]
    
    return (rho_s,r_s,a,b,g)
    
# DEKEL+ PROFILE
def get_params_dekel(rvir,mstar,mvir,fit_func=='new'):
    logMM=log10(mstar/mvir)
    
    if fit_func=='old':
        popt_s1=[1.76438231e-01, 1.71787327e-02, 1.42328182e+02, 1.23488025e+00]
        popt_c2=[1.86247747e+02, 1.37125654e+00, 0.00000000e+00, 1.41574162e-01]

        c2_DMO=10**(1.025-0.097*np.log10(mvir*0.671/1e12))
        s1_DMO=(1.+0.03*c2_DMO)/(1.+0.01*c2_DMO)

        s1_ratio=s1_function1(mstar/mvir,*popt_s1)
        s1=s1_ratio*s1_DMO

        c2_ratio=exp1min_func(mstar/mvir,*popt_c2)
        c2=c2_ratio*c2_DMO

        a=(1.5*s1-2.*(3.5-s1)*sqrt(0.01)*sqrt(c2))/(1.5-(3.5-s1)*sqrt(0.01)*sqrt(c2))
        c=((s1-2.)/((3.5-s1)*sqrt(0.01)-1.5/sqrt(c2)))**2
        
    if fit_func=='new':
        popt_s1=[1.30278262e-03,2.86267536e+00,3.20335633e-01]
        popt_c2=[1.86247747e+02, 1.37125654e+00, 0.00000000e+00, 1.41574162e-01]
        
        c2_DMO=10**(1.025-0.097*log10(mvir*0.671/1e12))
        s1_DMO=(1+0.03*c2_DMO)/(1+0.01*c2_DMO)
        
        s1_ratio=s1_zhao_dmo(mstar/mvir,*popt_s1)
        s1=s1_ratio*s1_DMO
        
        c2_ratio=exp1min_func(mstar/mvir,*popt_c2)
        c2=c2_ratio*c2_DMO
        
        a=(1.5*s1-2.*(3.5-s1)*sqrt(0.01)*sqrt(c2))/(1.5-(3.5-s1)*sqrt(0.01)*sqrt(c2))
        c=((s1-2.)/((3.5-s1)*sqrt(0.01)-1.5/sqrt(c2)))**2
        
    return (c, a, 2, 3, rvir, mvir)    

def s1_zhao_dmo(x,x0,nu,spp):
    return 1./(1.+pow(x/x0,nu))+spp*log10(1.+pow(x/x0,nu))

def s1_function1(x,x1,x2,nu1,nu2):#,eta):
    return 1.+np.log10(pow(1.+x/x1,-nu1)+pow(x/x2,nu2))

def exp1min_func(x,c1,nu,delta,mu):
    return 1.+c1*pow(x,nu)-pow(x,mu)-c1*pow(1.e-6,nu)+pow(1.e-6,mu)
    
###################################################################
# CONVERSION BETWEEN (a,c), (s1,c2), and (s1,cmax) for s1=s(r1)

# CONCENTRATION PARAMETER CORRESPONDING TO s(r2)=2 []
def c2_function(params, model='dekel'):

    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params
        return c*(1.5/(2.-a))**2

    if model == 'einasto':
        (c, n, Rvir, Mvir) = params   
        return c
        
    if model == 'enfw':
        (c,a,Rvir,Mvir) = params
        return c/(2.-a)   

def params_from_s1cmax(s1,cmax,r1,rvir,mvir):
    x12=np.sqrt(r1/rvir)
    c12=np.sqrt(cmax)
    a=(s1-2.*(3.5-s1)*x12*c12)/(1.-(3.5-s1)*x12*c12)
    c=((s1-2.)/((3.5-s1)*x12*c12-1./c12))**2
    return (c, a, 2, 3, rvir, mvir)

def params_from_s1c2(s1,c2,r1,rvir,mvir):
    x12=np.sqrt(r1/rvir)
    c12=np.sqrt(c2)
    a=(1.5*s1-2.*(3.5-s1)*x12*c12)/(1.5-(3.5-s1)*x12*c12)
    c=((s1-2.)/((3.5-s1)*x12-1.5/c12))**2.
    return (c, a, 2, 3, rvir, mvir)

###################################################################
# LENSING PROPERTIES (NUMERICAL INTEGRATIONS)

# PROJECTED SURFACE DENSITY [Msun kpc^-2]
def surfdens(r, params,truncated=True):
    (c, a, b, g, Rvir, Mvir) = params
    rc = float(Rvir)/c
    x = r/rc
    rhoc=q_bar(c, a, b, g, Mvir, Rvir)*float(3.-a)/3.
    
    if truncated:
        xmax=c
    else:
        xmax=np.inf
    
    if np.size(r)==1:
        if x==c:
            return 0.
        else: 
            return 2.*rhoc*rc*quad(integrand_surfdens,x,xmax,args=(x,a))[0]
    else: 
        Sigma_array=np.nan*np.ones_like(r)
        for i in range(np.size(r)):
            if x[i]==c:
                Sigma_array[i]=0. 
            else:
                Sigma_array[i]=2.*rhoc*rc*quad(integrand_surfdens,x[i],xmax,args=(x[i],a))[0]
        return Sigma_array

def integrand_surfdens(x,X,a):
    return pow(x,1.-a)*pow(1.+pow(x,0.5),2.*a-7)/np.sqrt(pow(x,2)-pow(X,2))

# PROJECTED CUMULATIVE MASS [Msun]
def cummass(r,params,truncated=True):
    cummass=np.nan*np.ones_like(r)
    for i in range(np.size(r)):
        cummass[i]=2.*np.pi*quad(integrand_cummass,0.,r[i],args=(params,truncated))[0]
    return cummass

def integrand_cummass(r,params,truncated=True):
    return r*surfdens(r, params, truncated)
    
# CONVERGENCE [in virial units, i.e., Mvir/pi Rvir^2, Sigma_crit=1]
def convergence(r,params,truncated=True):
    (c, a, b, g, Rvir, Mvir) = params
    return cummass(r,params,truncated)/Mvir*pow(Rvir/r,2)

# SHEAR [in virial units, i.e., Mvir/pi Rvir^2, Sigma_crit=1]
def shear(r,params,truncated=True):
    (c, a, b, g, Rvir, Mvir) = params
    surfvir=Mvir/(np.pi*Rvir**2)
    sigmab=cummass(r,params,truncated)/Mvir*pow(Rvir/r,2)
    sigma=surfdens(r, params,truncated)/surfvir
    return sigmab-sigma
    
###################################################################
# AUXILIARY FUNCTIONS

def q_bar(c, a, b, g, Mvir, Rvir):
    brho_vir = Mvir/(4*np.pi/3.*pow(Rvir,3))
    return brho_vir*pow(c,a)*pow((1.+pow(c,1./b)),b*(g-a))

def mu(c, a, b, g):
    return pow(1.+pow(c,1./b),b*(g-a))/float(pow(c,3.-a))

def calc_total_mass(r, rho, Rmax):
    return sum(rho[r<=Rmax]*4*np.pi*r[r<=Rmax]**2*np.concatenate(([r[0]], np.diff(r[r<=Rmax]))))

def dM_dicintio(r,rho_s,r_s,a,b,g):
    return 4.*np.pi*rho_s*pow(r,2)/(pow(r/r_s,g)*pow(1.+pow(r/r_s,a),(b-g)/a))

def dM_gNFW(r,rho_s,r_s,g):
    return 4.*np.pi*rho_s*pow(r,2)/( pow(r/r_s,g) * pow(1.+r/r_s,3-g) )

def dM_einasto(r,q,h,n):
    return 4.*np.pi*pow(r,2) * q*np.exp( -pow( (r/h), 1./n ) )

def Binc(a,b,x):
    return betainc(a,b,x)*gamma(a)*gamma(b)/gamma(a+b)
    
def Bintegrand(t,a,b):
    return pow(t,a-1.)*pow(1-t,b-1.)

def Dbetainc(a,b,x1,x2=1):
    if np.size(x1)==1:
        return quad(Bintegrand,x1,x2,args=(a,b,))[0]
    elif np.size(x1)>1:
        Dbeta=np.nan*np.zeros_like(x1)
        for i in range(np.size(x1)):
            Dbeta[i]=quad(Bintegrand,x1[i],x2,args=(a,b,))[0]
        return Dbeta

def chi(x):
    u = pow(x,0.5)
    return u/(1.+u)

def B9(a,x):
    """
    Incomplete beta function B(a,b,x) with b=9
    betainc(a,9,x)
    """
    B=0.
    for i in range(9):
        B+=factorial(8)/factorial(i)*gamma(a)/gamma(a+9-i)*x**(a+8-i)*(1.-x)**i
    return B

###################################################################
# AUXILIARY FUNCTIONS FOR THE DISTRIBUTION FUNCTION 

def distribution_integrand(x,E,params,model='dekel'):
    (c, a, b, g, Rvir, Mvir) = params
    integrand= np.sqrt(E+U_adim(x, params, model='dekel'))*lambda_(x,params)
    return np.nan_to_num(integrand)

def U_adim(x, params, model='dekel'): # U in [GMvir/Rvir]

    if model == 'dekel': 
        (c, a, b, g, Rvir, Mvir) = params

        if (b==1 and g==3): 
            return c*mu(c,a,b,g)/float(2.-a)*(pow(x/np.array(1.+x),2-a)-pow(c/float(1.+c),2.-a))-1.
        
        if (b==2 and g==3):
            chi = pow(x,0.5)/np.array(1.+pow(x,0.5))
            chic = pow(c,0.5)/np.array(1.+pow(c,0.5))
            return -1.-2*c*mu(c,a,b,g)*((pow(chic,2*(2-a))-pow(chi,2*(2-a)))/float(2*(2-a))-(pow(chic,2*(2-a)+1)-pow(chi,2*(2-a)+1))/float(2*(2-a)+1))
        else:
            dU = lambda y: M_function(y,params,model)/Mvir/array(pow(y,2))
            return array([-1.-quad(dU,xt,c)[0] for xt in x])

def lambda_(x,params):
    return (g_prim(x,params)*dpsi(x,params)-dnudpsi(x,params)*d2psi(x,params))/pow(dpsi(x,params),2)

def g_prim(x,params):
    (c, a, b, g, Rvir, Mvir) = params
    return (-3+a)*c**2*(48*a+3*(35+44*a)*x**0.5+4*x*(77+24*a)+245*x**1.5)/32./np.pi/(x+x**1.5)**4

def dnudpsi(x,params):
    (c, a, b, g, Rvir, Mvir) = params
    factor=(3.-a)*pow(c,2)/(4.*np.pi)
    numerator=2.*a+(5.25+3.*a)*pow(x,0.5)+8.75*x
    denominator=x**3*pow(1.+pow(x,0.5),3)
    return factor*numerator/denominator

def d2psi(x,params):
    (c, a, b, g, Rvir, Mvir) = params
    return c*mu(c, a, b, g)*pow(1.+1./pow(x,0.5),2*a)*(-1+a+2*pow(x,0.5))/pow(1.+pow(x,0.5),7)

def dpsi(x,params):
    (c, a, b, g, Rvir, Mvir) = params
    return -c*mu(c, a, b, g)/pow(x,a-1)/pow(1.+pow(x,0.5),6-2*a)

###################################################################



