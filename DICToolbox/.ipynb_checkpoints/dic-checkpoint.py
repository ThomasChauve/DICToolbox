# -*- coding: utf-8 -*-
'''
Created on 21 août 2015

.. py:module:: dic

Toolbox for data obtained using 7D software create by Pierre Vacher
7D software is software to perform Digital Image Correlation (DIC). This toolbox is used to analysed the data export from 7D software by using the exportation by 'txt>>gdr'

@author: Thomas Chauve
@contact: thomas.chauve@lgge.obs.ujf-grenoble.fr
@license: CC-BY-CC
'''

import DICToolbox.image2d as im2d
import DICToolbox.symetricTensorMap as sTM

import numpy as np
from skimage import io
import os
import pylab
import matplotlib.pyplot as plt

import matplotlib.cm as cm

class dic(object):
    '''
    'dic' is a python class to analyse data from 2D digital image correlation
    
    .. note:: this toolbox is build to work with the output from 7D using txt>>gdr export   
    '''
    pass


    def __init__(self, time, u, v, strain, oxy, micro, grains):
        '''
        :param adr_data: folder path where the 7D output are store
        :type adr_data: str
        :param resolution: pixel size of the image used for DIC (millimeters)
        :type resolution: float
        :param time_step: time step between the picture (seconds)
        :type time_step: float
        :param adr_micro: path for the black and white skeleton of the microstructure (bmp format) (default 0 - no microstructure)
        :type adr_micro: str
        '''
        
        # include time_step in the object
        self.time=time
        self.u=u
        self.v=v
        self.strain=strain
        self.oxy=oxy
        self.micro=micro
        self.grains=grains
        
    def crop(self,xmin=0,xmax=0,ymin=0,ymax=0):
        '''
        Crop function to select the area of interest
        
        :return: crop dic object
        :rtype: dic
        :Exemple: >>> data.crop()
        
        .. note:: clic on the top left corner and bottom right corner to select the area
        '''
        
        if (xmin+xmax+ymin+ymax)==0:
            # plot the data
            h=self.strain[len(self.u)-1].t22.plot()
            self.micro.plotBoundary()
            # select top left and bottom right corner for crop
            print('Select top left and bottom right corner for crop :')
            x=np.array(pylab.ginput(2))/self.u[1].res
            plt.close("all")
            # create x and Y coordinate

            xx=[x[0][0],x[1][0]]
            yy=[x[0][1],x[1][1]]
            # size of the initial map
            ss=np.shape(self.u[0].field)
            # find xmin xmax ymin and ymax
            xmin=np.int(np.ceil(np.min(xx)))
            xmax=np.int(np.floor(np.max(xx)))
            ymin=np.int(ss[0]-np.ceil(np.max(yy)))
            ymax=np.int(ss[0]-np.floor(np.min(yy)))

        
        for i in list(range(len(self.strain))):
            # crop the map
            self.u[i].field=self.u[i].field[ymin:ymax, xmin:xmax]
            self.v[i].field=self.v[i].field[ymin:ymax, xmin:xmax]
            #
            ss=np.shape(self.u[i].field)
            #self.exx[i].field=self.exx[i].field[ymin:ymax, xmin:xmax]
            #self.eyy[i].field=self.eyy[i].field[ymin:ymax, xmin:xmax]
            #self.exy[i].field=self.exy[i].field[ymin:ymax, xmin:xmax]
            #self.eeq[i].field=self.eeq[i].field[ymin:ymax, xmin:xmax]
            self.strain[i].t11.field=self.strain[i].t11.field[ymin:ymax, xmin:xmax]
            self.strain[i].t22.field=self.strain[i].t22.field[ymin:ymax, xmin:xmax]
            #self.strain[i].t33.field=self.strain[i].t33.field[ymin:ymax, xmin:xmax]
            self.strain[i].t12.field=self.strain[i].t12.field[ymin:ymax, xmin:xmax]
            #self.strain[i].t13.field=self.strain[i].t13.field[ymin:ymax, xmin:xmax]
            #self.strain[i].t23.field=self.strain[i].t23.field[ymin:ymax, xmin:xmax]
            self.strain[i].t33.field=np.zeros(ss)
            self.strain[i].t13.field=np.zeros(ss)
            self.strain[i].t23.field=np.zeros(ss)
            
            self.oxy[i].field=self.oxy[i].field[ymin:ymax, xmin:xmax]
            
        self.micro.field=self.micro.field[ymin:ymax, xmin:xmax]
        self.grains=self.micro.grain_label()
        
        # replace grains boundary with NaN number
        self.grains.field=np.array(self.grains.field,float)
        idx=np.where(self.micro.field==1)
        self.grains.field[idx]=np.nan
        
        return np.array([xmin,xmax,ymin,ymax])
        
    def strain_macro(self,nb_line=3):
        '''
        Compute the macroscopic strain along the y-axis (vertical) using 2 methods :
        
        (1)-mean of eyy \n
        (2)-dic line (see Fanny Grennerat Thesis)
            
        :param nb_line: number of line skip at the top and bottom part of the picture during (2) evalutation
        :type nb_line: int
        :return: macro_eyy, method (1)
        :rtype: array
        :return: macro_line, method (2)
        :rtype: array
        :return: time, time of each point (seconds)
        :rtype: array
        :Exemple:
            >>> [time,macro_eyy,macro_line]=data.strain_macro()
            >>> plt.plot(time,macro_eyy)
            >>> plt.show()
            
        .. warning:: If you have macro_eyy~=-macro_line, it is because you did not correct the gdr output 
        '''
        
        time=np.zeros(len(self.u))
        macro_eyy=np.zeros(len(self.u))
        macro_line=np.zeros(len(self.u))
        # number of line along y
        nby=np.shape(self.v[0].field)[0]-1
        for i in list(range(len(self.u))):
            time[i]=(i+1)*self.time
            macro_eyy[i]=np.nanmean(self.strain[i].t22.field)
            macro_line[i]=(np.nanmean(self.v[i].field[nb_line,:])-np.nanmean(self.v[i].field[nby-nb_line,:]))/((nby-2*nb_line)*self.v[i].res)
            
            
        return time,macro_eyy,-macro_line
    
    def find_pic(self,strain_step=0.001,minst=0):
        '''
        Find picture number as between 2 picture the increment of deformation is macro_strain
        
        :param strain_step: macroscopic strain between 2 picture.
        :type strain_step: float
        :return: nb_img
        :rtype: array
        '''
        
        # extract macroscopic strain from the data
        [time,macro_eyy,macro_line]=self.strain_macro()
        
        nb_img=[]
        
        id=np.isnan(macro_line)
        macro_line[id]=0
        
        #last strain record
        if minst==0:
            minstrain=macro_line[len(macro_line)-1]
        else:
            minstrain=minst
        
        
        #number of iteration
        it=np.int32(minstrain/strain_step)

        for i in list(range(-it)):
            # find where the macro strain is higher that the step
            idx=np.where(macro_line>-1.*np.float(i+1)*strain_step)
            # take the minimum
            label=np.where(macro_line==np.min(macro_line[idx]))
            nb_img.append(list(label[0]))
            
            
        res=np.array(nb_img)    
        #np.savetxt('wanted_image.txt',res,fmt='%i')
            
        return res
    
    def noise_map(self):
        '''
        compute the variance of a set of map for u, v, exx, eyy, exy, eeq
        :return su,sv,sexx,seyy,sexy,seeq:
        :rtype su,sv,sexx,seyy,sexy,seeq: im2d.image2d
        '''
        l, c = np.shape(self.u[0].field)
        st11 = np.zeros((l, c, len(self.u)))
        st22 = np.zeros((l, c, len(self.u)))
        st12 = np.zeros((l, c, len(self.u)))
        steq = np.zeros((l, c, len(self.u)))
        su = np.zeros((l, c, len(self.u)))
        sv = np.zeros((l, c, len(self.u)))
        
        for i in list(range(len(self.u))):
            st11[:,:,i]=self.strain[i].t11.field
            st22[:,:,i]=self.strain[i].t22.field
            st12[:,:,i]=self.strain[i].t12.field
            su[:,:,i]=self.u[i].field
            sv[:,:,i]=self.v[i].field
            steq[:,:,i]=self.strain[i].eqVonMises().field
            
        return im2d.image2d(np.nanstd(su,axis=2),self.micro.res), im2d.image2d(np.nanstd(sv,axis=2),self.micro.res), im2d.image2d(np.nanstd(st11,axis=2),self.micro.res), im2d.image2d(np.nanstd(st22,axis=2),self.micro.res), im2d.image2d(np.nanstd(st12,axis=2),self.micro.res), im2d.image2d(np.nanstd(steq,axis=2),self.micro.res)
        
