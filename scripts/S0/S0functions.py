import numpy as np

def fit_OrnsteinZernike(x, s0, xi):
    # S(k) = S(0)/(1+xi*k^2)
    return s0/(1.+xi*x**2.)

def ideal_mu(c, c0, kbt):
    return kbt*np.log(c/c0)

def get_activity_coefficient(saa,sab,nratio):
    return 1./(saa-sab*(nratio)**0.5)

def excess_mu(ca,saa,sab,nratio,c0,kbT):
    # ca, nratio are 1-D array with N-dim
    # saa,sab are (2,N) array, [S, S_error]
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    
    integrand = 1./(saa[:,0]-sab[:,0]*(nratio)**0.5) -1.
    integrand_error = 0.5*abs(1./((saa[:,0]-saa[:,1])-(sab[:,0]+sab[:,1])*(nratio)**0.5)-1./((saa[:,0]+saa[:,1])-(sab[:,0]-sab[:,1])*(nratio)**0.5))
    logca = np.log(ca)
    
    mask = ~np.isnan(integrand)
    int_func = interp1d(logca[mask], integrand[mask], kind='linear',fill_value='extrapolate')
    int_error_func = interp1d(logca[mask], integrand_error[mask], kind='linear',fill_value='extrapolate')
    
    # integrate
    ca_mu = []
    n_sample = 0
    for i,a in enumerate(logca):
        if ca[i] > 0:
            n_sample+=1
            mu_now = kbT*quad(int_func, np.log(c0), a)[0]
            if n_sample >= 2:
                mu_error_now = kbT*quad(int_error_func, np.log(c0), a)[0]/(n_sample-1)**0.5
            else:
                mu_error_now = kbT*quad(int_error_func, np.log(c0), a)[0]
            ca_mu.append([ca[i],mu_now, mu_error_now])
        
    return np.asarray(ca_mu)

def get_betadmua_dlnxa(saa,sab,sbb,xa):
    xb = 1. - xa
    return 1./(xb*saa+xa*sbb-2*(xa*xb)**0.5*sab)

def excess_mu_x_GH(saa,sab,sbb,xa,x0,kbT):
    # xa are 1-D array with N-dim
    # saa,sab,sbb are (2,N) array, [S, S_error]
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    
    integrand = get_betadmua_dlnxa(saa[:,0],sab[:,0],sbb[:,0],xa)
    integrand_error = np.abs(get_betadmua_dlnxa(saa[:,0],sab[:,0],sbb[:,0],xa) \
                             - get_betadmua_dlnxa(saa[:,0]+saa[:,1],sab[:,0]-sab[:,1],sbb[:,0]+sbb[:,1],xa))
    logxa = np.log(xa)
    
    mask = ~np.isnan(integrand)
    int_func = interp1d(logxa[mask], integrand[mask], kind='linear',fill_value='extrapolate')
    int_error_func = interp1d(logxa[mask], integrand_error[mask], kind='linear',fill_value='extrapolate')
    
    # integrate
    xa_mu = []
    n_sample = 0
    for i,a in enumerate(logxa):
        n_sample+=1
        mu_now = kbT*quad(int_func, np.log(x0), a)[0]
        if n_sample >= 2:
            mu_error_now = kbT*quad(int_error_func, np.log(x0), a)[0]/(n_sample-1)**0.5
        else:
            mu_error_now = kbT*quad(int_error_func, np.log(x0), a)[0]
        xa_mu.append([xa[i],mu_now, mu_error_now])
        
    return np.asarray(xa_mu)
