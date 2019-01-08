#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
.. py:module:: correct7Dgdr

.. warning:: used this program on the data export using 'txt>>gdr' of the 7D version 0.9.5.178

This program correct the error in the output data from 7D version 0.9.5.178 were dy=-dy and exy=-exy

:Example: 
    >> run correct7Dgdr in the folder contening the output '.txt'
'''

import numpy as np
import os


def correction(adr)
# find the curent directory
current_dir = os.getcwd()

folder='CorrectedData/'
os.mkdir(folder)

for file in os.listdir(current_dir):
    if file.endswith(".txt"):
        #print(file)
        # load the file
        data=np.loadtxt(file)
        # correct the data
        data[:,3]=-data[:,3]
        data[:,7]=-data[:,7]
        # save the file
        np.savetxt(folder+file, data)