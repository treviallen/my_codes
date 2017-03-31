from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack



def quadratic_vertex(x, y, vertex, *kwargs):
    
    def fit_quadratic_vertex(c, x):
    
        xx = x - vertex
        
        return c[0] * xx**2 + c[1]
    
    c0*(x-vertex)**2 + c1