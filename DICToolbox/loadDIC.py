import DICToolbox.dic as dic
import DICToolbox.image2d as im2d
import DICToolbox.symetricTensorMap as sTM
from skimage import io

import os
import numpy as np

def load7D(adr_data, resolution, time_step, adr_micro=0):
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
        time=time_step
                
        # find out file from 7D
        output=os.listdir(adr_data)
        output.sort()
        # loop on all the output
        u=[]
        v=[]
        strain=[]
        oxy=[]
        for i in list(range(len(output))):
            # load the file
            data=np.loadtxt(adr_data + '/' + output[i])
            # find not indexed point
            id=np.where(data[:,4]!=0.0)
            # replace not indexed point by NaN value
            data[id,2:4]=np.NaN
            data[id,5:8]=np.NaN
            # for the first step it extract size sx,sy and window correlation value used in 7D
            if (i==0):
                n_dic=np.abs(data[0,1]-data[1,1])
                nb_pix=np.size(data[:,0])
                sx=np.int32(np.abs(data[0,0]-data[nb_pix-1,0])/n_dic+1)
                sy=np.int32(np.abs(data[0,1]-data[nb_pix-1,1])/n_dic+1)
            # Build image at time t=i
            u.append(im2d.image2d(np.transpose(np.reshape(data[:,2]*resolution,[sx,sy])),n_dic*resolution))
            v.append(im2d.image2d(np.transpose(np.reshape(data[:,3]*resolution,[sx,sy])),n_dic*resolution))
            exx=im2d.image2d(np.transpose(np.reshape(data[:,5],[sx,sy])),n_dic*resolution)
            eyy=im2d.image2d(np.transpose(np.reshape(data[:,6],[sx,sy])),n_dic*resolution)
            exy=im2d.image2d(np.transpose(np.reshape(data[:,7],[sx,sy])),n_dic*resolution)
            if i==0:
                eaz=im2d.image2d(np.zeros([sy,sx]),n_dic*resolution)
                
            strain.append(sTM.symetricTensorMap(exx,eyy,eaz,exy,eaz,eaz))
            oxy.append((u[i].diff('y')-v[i].diff('x'))*.5)
        
        minx=np.int32(np.min(data[:,0])/n_dic)-1
        maxx=np.int32(np.max(data[:,0])/n_dic)
        miny=np.int32(np.min(-data[:,1])/n_dic)-1
        maxy=np.int32(np.max(-data[:,1])/n_dic)
                
        # open the microstructure
        if adr_micro==0:
            micro=im2d.micro2d(np.zeros([sy,sx]),n_dic*resolution)
        else:
            micro_bmp = io.imread(adr_micro)
            micro=im2d.micro2d(micro_bmp[miny:maxy,minx:maxx,0]/np.max(micro_bmp),n_dic*resolution)
        
        grains=micro.grain_label()
        
        # replace grains boundary with NaN number
        grains.field=np.array(grains.field,float)
        idx=np.where(micro.field==1)
        grains.field[idx]=np.nan
        
        return dic.dic(time, u, v, strain, oxy, micro, grains)
    
    
    
def loadDICe(adr_data, resolution, time_step, adr_micro=0,strain_com=1):
        '''
        :param adr_data: folder path where the 7D output are store
        :type adr_data: str
        :param resolution: pixel size of the image used for DIC (millimeters)
        :type resolution: float
        :param time_step: time step between the picture (seconds)
        :type time_step: float
        :param adr_micro: path for the black and white skeleton of the microstructure (bmp format) (default 0 - no microstructure)
        :type adr_micro: str
        :param strain_com: Computation method for strain calculation 0-small deformation 1-Green Lagrange (default 1)
        :type strain_com: bool
        '''
        
        # include time_step in the object
        time=time_step
                
        # find out file from 7D
        output=os.listdir(adr_data+'results/')
        output.sort()
        output=output[2:-1]
        # loop on all the output
        u=[]
        v=[]
        strain=[]
        oxy=[]
        for i in list(range(len(output))):
            # load the file
            data=np.loadtxt(adr_data+'results/'+ output[i],skiprows=1,delimiter=',',usecols=[1,2,3,4,9])
            # find not indexed point
            id=np.where(data[:,-1]==0)
            # replace not indexed point by NaN value
            data[id[0],2]=np.NaN
            data[id[0],3]=np.NaN
            # for the first step it extract size sx,sy and window correlation value used in 7D
            if (i==0):
                n_dic=np.abs(data[0,1]-data[1,1])
                nb_pix=np.size(data[:,0])
                sx=np.int32(1+np.abs(np.min(data[:,0])-np.max(data[:,0]))/n_dic)
                sy=np.int32(1+np.abs(np.min(data[:,1])-np.max(data[:,1]))/n_dic)
            # Build image at time t=i
            tmp_u=np.zeros([sx,sy])
            tmp_u[:,:]=np.NaN
            tmp_v=np.zeros([sx,sy])
            tmp_v[:,:]=np.NaN
            # find x y position in the table
            x=(data[:,0]-np.min(data[:,0]))/n_dic
            y=(data[:,1]-np.min(data[:,1]))/n_dic
            # include data
            tmp_u[x.astype(int),y.astype(int)]=data[:,2]
            tmp_v[x.astype(int),y.astype(int)]=data[:,3]
            # Create im2d data object u and v
            imu=im2d.image2d(np.transpose(tmp_u*resolution),n_dic*resolution)
            imv=im2d.image2d(np.transpose(tmp_v*resolution),n_dic*resolution)
            # compute Strain
            if strain_com==0:
                exx=imu.diff('x')
                eyy=imv.diff('y')
                exy=(imu.diff('y')+imv.diff('x'))*0.5
                tmp_oxy=(imu.diff('y')-imv.diff('x'))*0.5
            elif strain_com==1:
                exx=(imu.diff('x')+(imu.diff('x').pow(2)+ imv.diff('x').pow(2))*0.5)
                eyy=(imv.diff('y')+(imu.diff('y').pow(2)+ imv.diff('y').pow(2))*0.5)
                exy=((imu.diff('y')+imv.diff('x'))*0.5+(imu.diff('x')*imu.diff('y')+imv.diff('x')*imv.diff('y'))*0.5)
                tmp_oxy=(imu.diff('y')-imv.diff('x'))*0.5#+(imu.diff('x')*imu.diff('y')-imv.diff('x')*imv.diff('y'))*0.5

    
            u.append(imu)
            v.append(imv)
            # faxe componante to use 3D strain tensor object
            if i==0:
                eaz=im2d.image2d(np.zeros([sy-1,sx-1]),n_dic*resolution)
            
            exx.field=-exx.field
            eyy.field=-eyy.field
            exy.field=-exy.field
            strain.append(sTM.symetricTensorMap(exx,eyy,eaz,exy,eaz,eaz))
            oxy.append(tmp_oxy)
        
        minx=np.int32(np.min(data[:,0])/n_dic)-1
        maxx=np.int32(np.max(data[:,0])/n_dic)
        miny=np.int32(np.min(-data[:,1])/n_dic)-1
        maxy=np.int32(np.max(-data[:,1])/n_dic)
        
        # open the microstructure
        if adr_micro==0:
            micro=im2d.micro2d(np.zeros([sy,sx]),n_dic*resolution)
        else:
            micro_bmp = io.imread(adr_micro)
            micro=im2d.micro2d(micro_bmp[miny:maxy,minx:maxx,0]/np.max(micro_bmp),n_dic*resolution)
        
        grains=micro.grain_label()
        
        # replace grains boundary with NaN number
        grains.field=np.array(grains.field,float)
        idx=np.where(micro.field==1)
        grains.field[idx]=np.nan
        
        return dic.dic(time, u, v, strain, oxy, micro, grains)