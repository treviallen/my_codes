from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack



def quadratic_vertex(x, y, vertex):
    from scipy.odr import Data, Model, ODR, models
    import scipy.odr.odrpack as odrpack

    def fit_quadratic_vertex(c, x):
    
        xx = x - vertex
        
        return c[0] * xx**2 + c[1]
    
    c0*(x-vertex)**2 + c1
    
#from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x >= hx
    xmod[idx] = 1
    return xmod
    
def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    #idx1 = x <= hx
    #idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)
    
    return ans1 + ans2
    
def get_linear_fixed_slope(slope, x, y):
    def linear_fixed_slope(c, x):  
                                       
        return meanslope * x + c[0]
    
    # do something here
    
        
    return intercept


# function for linear ODR
def lin_odr(c,x):
    return c[0] * x + c[1]

def odr_lin_reg(xin, yin, start):
    '''
    start is starting conditions, i.e. [m, c]
    xin, yin = input arrays    
    '''
    data = odrpack.RealData(xin, yin)
    
    lin_reg = odrpack.Model(lin_odr)
    odr = odrpack.ODR(data, lin_reg, beta0=start)
    odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
    out = odr.run()
    out.pprint()
    
    # set coeffs
    m = out.beta[0]
    c = out.beta[1]

    return m, c

def lsq_lin_reg(xin, yin, start):
    '''
    start is starting conditions, i.e. [m, c]
    xin, yin = input arrays    
    '''
    data = odrpack.RealData(xin, yin)
    
    lin_reg = odrpack.Model(lin_odr)
    odr = odrpack.ODR(data, lin_reg, beta0=start)
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
    out = odr.run()
    out.pprint()
    
    # set coeffs
    m = out.beta[0]
    c = out.beta[1]

    return m, c

# From: https://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb
def lsq_reg_with_confidence(x, y):
    """Linear fit of x and y with uncertainty and plots results."""
    
    import numpy as np
    import scipy.stats as stats
    
    x, y = np.asarray(x), np.asarray(y)
    n = y.size
    '''
    yerr = np.abs(4*np.random.randn(len(x))) + 2
    p, cov = np.polyfit(x, y, 1, w=1/yerr, cov=True)  # coefficients and covariance matrix
                             # evaluate the polynomial at x
    perr = np.sqrt(np.diag(cov))     # standard-deviation estimates for each coefficient
    R2 = np.corrcoef(x, y)[0, 1]**2  # coefficient of determination between x and y
    
    chi2red = np.sum((resid/yerr)**2)/(n - 2)  # Chi-square reduced
    '''
    p = np.polyfit(x, y, 1)  # coefficients and covariance matrix
    yfit = np.polyval(p, x)
    resid = y - yfit
    s_err = np.sqrt(np.sum(resid**2)/(n - 2))  # standard deviation of the error (residuals)
    
    # Confidence interval for the linear fit:
    t = stats.t.ppf(0.975, n - 2)
    ci = t * s_err * np.sqrt(1/n + (x - np.mean(x))**2/np.sum((x-np.mean(x))**2))
    # Prediction interval for the linear fit:
    pi = t * s_err * np.sqrt(1 + 1/n + (x - np.mean(x))**2/np.sum((x-np.mean(x))**2))

    '''
    # Plot
    plt.fill_between(x, yfit+pi, yfit-pi, color=[1, 0, 0, 0.1], edgecolor='')
    plt.fill_between(x, yfit+ci, yfit-ci, color=[1, 0, 0, 0.15], edgecolor='')
    plt.errorbar(x, y, yerr=yerr, fmt = 'bo', ecolor='b', capsize=0)
    plt.plot(x, yfit, 'r', linewidth=3, color=[1, 0, 0, .8])
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.title('$y = %.2f \pm %.2f + (%.2f \pm %.2f)x \; [R^2=%.2f,\, \chi^2_{red}=%.1f]$'
              %(p[1], perr[1], p[0], perr[0], R2, chi2red), fontsize=20, color=[0, 0, 0])  
    plt.xlim((0, n+1))
    plt.show()
    '''
    # return confidence intervals
    return yfit+ci, yfit-ci
"""    
##############################################################
# example
##############################################################
'''
x & y are some data we want to regress!

xrng & yrng are some target model

sx & sy are arrays of sigma values corresponding to x & y

'''

data = odrpack.RealData(x, y)
# or
data = odrpack.RealData(x, y, sx=sx, sy=sy) # sx & sy are optional!


bilin_reg = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data, bilin_reg, beta0=[1.0, -3.0, 1.0, 8.6])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
out.pprint()

a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3] # x hinge point

yrng = b + a * xrng
yhinge = b + a * hx
idx = xrng > hx
yrng[idx] = c * (xrng[idx]-hx) + yhinge
"""