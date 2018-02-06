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