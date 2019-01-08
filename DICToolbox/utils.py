import numpy as np
import os

def fmt(x, pos):
    '''
    function to format tick in the figure as x*10^y
    '''
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
