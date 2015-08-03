#!/usr/bin/env python

"""Newton-Raphson solver using fixed-points and numerical first-derivative
approximation.

"""

def iterative_improve(good_enough, improve):

    """ """
    def iterate(guess):
        next_try = improve(guess)
        if good_enough(guess, next_try):
            return next_try
        else:
            return iterate(next_try)
    return lambda x: iterate(x)

def fixed_point(f, first_guess, epsilon=1.0e-12):
    """ """
    #def good_enough(v1, v2):
    #    return abs(v1-v2)/v2 < epsilon
    return iterative_improve(good_enough, f)(first_guess)

def good_enough(v1, v2, epsilon=1e-8):
    return abs( (v1 - v2) / v2 ) < epsilon
        

def approx_deriv(f, dx=1.0e-8):
    """
    
    """
    return lambda x: float(f(x + dx) - f(x)) / dx

def newton_transform(f):
    """"""
    return lambda x: x - float(f(x)) / ( approx_deriv(f) (x) )
    
    
def newton_raphson(f, guess):
    """
    newton_raphson(f, guess)

    Newton Raphson procedure
    
    Parameters
    ----------
    f : function
        
    Returns
    -------
    float
        Zero of the specified function.
    """
    return fixed_point(newton_transform(f), guess)


#if __name__ == '__main__':
#    main()
