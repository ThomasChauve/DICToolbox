import numpy as np
import os

def correction7Ddata(adr,outputfolder='CorrectedData/'):
               
    os.mkdir(outputfolder)

    for file in os.listdir(adr):
        if file.endswith(".txt"):
            # load the file
            data=np.loadtxt(adr+file)
            # correct the data
            data[:,3]=-data[:,3]
            data[:,7]=-data[:,7]
            # save the file
            np.savetxt(outputfolder+file, data)

def fmt(x, pos):
    '''
    function to format tick in the figure as x*10^y
    '''
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)