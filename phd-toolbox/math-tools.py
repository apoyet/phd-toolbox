import numpy as np

# Linear Regression

def lin_reg(x,y):
    '''
    Simple linear regression. Takes as input 2 arrays x and y being thte coordinates of the 
    original data. 
    Return the coefficients a and b such as y = ax+b.
    
    Args:
        x: x coordinate of the original data
        y: y coordinate of the original data
    
    Returns: 
        a: float, slope of the linear regression
        b: float, origin of the linear regression
        [a,b] is returned as a list
    '''
    a = (np.sum((x - np.mean(x))*(y - np.mean(y))))/(np.sum((x - np.mean(x))**2))
    b = np.mean(y) - a*np.mean(x)
    return [a, b]
