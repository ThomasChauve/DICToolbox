# -*- coding: utf-8 -*-
'''
Created on 3 juil. 2015

image2d is a class used to manipulate image under matrix shape and to do the analyses on the picture

.. note:: It has been build to manipulate both aita data and dic data
.. warning:: As all function are applicable to aita and dic data, please be careful of the meaning of what you are doing depending of the input data used !  

@author: Thomas Chauve
@contact: thomas.chauve@univ-grenoble-alpes.fr
@license: CC-BY-CC
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from skimage import io
from skimage import morphology
import scipy
import pylab
import datetime
import math


class image2d(object):
    '''
    image2d is a class for map of scalar data
    '''
    pass

    def __init__(self, field, resolution):
        '''
        Constructor : field is a matrix and resolution is the step size in mm
        
        :param field: tabular of scalar data
        :type field: array
        :param resolution: step size resolution (millimeters) 
        :type resolution: float
        '''
        
        self.field=field
        self.res=resolution
        
    def plot(self,vmin=np.NaN,vmax=np.NaN,colorbarcenter=False,colorbar=cm.viridis):
        '''
        plot the image2d
        
        :param vmin: minimum value for the colorbar
        :type vmin: float
        :param vmax: maximun value for the colorbar
        :type vmax: float
        :param colorbarcenter: do you want center the colorbar around 0
        :type colorbarcenter: bool
        :param colorbar: colorbar from matplotlib.cm
        
        .. note:: colorbar : cm.jet for eqstrain-stress
        '''
        if np.isnan(vmin):
            vmin=np.nanmin(self.field)
            
        if np.isnan(vmax):
            vmax=np.nanmax(self.field)
        
        # size of the image2d
        ss=np.shape(self.field)
        # create image
        plt.imshow(self.field,aspect='equal',extent=(0,ss[1]*self.res,0,ss[0]*self.res),cmap=colorbar,vmin=vmin,vmax=vmax)
        
        if colorbarcenter:
            zcolor=np.max(np.max(np.abs(self.field)))
            plt.clim(-zcolor, zcolor)
        
        # set up colorbar
        plt.colorbar()
        
    def extract_data(self,pos=[]):
        '''
        Extract the value at the position 'pos' or where you clic
        
        :param pos: array [x,y] position of the data, if pos==[], clic to select the pixel
        :type pos: array
        '''
        
        if pos==[]:
            plt.imshow(self.field,aspect='equal')
            plt.waitforbuttonpress()
            print('select the pixel :')
            #grain wanted for the plot
            id=np.int32(np.array(pylab.ginput(1)))
        else:
            id=pos

        plt.close()    
        return self.field[id[0,1],id[0,0]],id
    
    def triple_junction(self):
        '''norm=matplotlib.colors.LogNorm()
        Localized the triple junction
        '''
        ss=np.shape(self.field)
        triple=[]
        pos=[]
        for i in list(range(ss[0]-2)):
            for j in list(range(ss[1]-2)):
                sub=self.field[i:i+2,j:j+2]
                id=np.where(sub[:]==sub[0,0])
                i1=len((id[0]))
                if (i1<(3)):
                    id=np.where(sub[:]==sub[0,1])
                    i2=len((id[0]))
                    if (i2<(3)):
                        id=np.where(sub[:]==sub[1,1])
                        i3=len((id[0]))
                        if ((i3==1 or i2==1 or i1==1) and i3<3):
                            triple.append(sub)
                            pos.append([i,j])
                            
        c=np.array(pos)
        z=np.arange(len(c[:,0]))                    
        plt.imshow(self.field)
        plt.plot(c[:,1],c[:,0],'+')
        
        for label, x, y in zip(z, c[:, 1], c[:, 0]):
            plt.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        return triple,c
    
    def imresize(self,res):
        '''
        Resize the image with nearest interpolation to have a pixel of the length given in res
        
        :param res: new resolution map wanted (millimeters)
        :type res: float
        '''
        # Fraction of current size
        zoom=float(self.res/res)
        self.res=res
        self.field=scipy.ndimage.interpolation.zoom(self.field,zoom,order=0,mode='nearest')
        
    
    def diff(self,axis):
        '''
        Derive the image along the axis define
        
        :param axis: 'x' to do du/dx or 'y' to do du/dy
        :type axis: str
        :return: differential of the map in respect to the 'axis' wanted
        :rtype: image2d
        
        .. warning:: 'x' axis is left to right and 'y' is bottom to top direction        
        '''
        if (axis=='x'):
            dfield=np.diff(self.field,axis=1)/self.res
            nfield=dfield[1:,:] # remove one line along y direction to have same size for both diff
        elif (axis=='y'):
            dfield=-np.diff(self.field,axis=0)/self.res # the - is from the convention y axis from bottom to top
            nfield=dfield[:,1:]
        else:
            print('axis not good')
            
        dmap=image2d(nfield,self.res)
        
        return dmap
    
    def __add__(self, other):
        '''
        Sum of 2 maps
        '''
        if (type(other) is image2d):
            return image2d(self.field+other.field,self.res)
        elif (type(other) is float):
            return image2d(self.field+other,self.res)
            
    def __sub__(self, other):
        '''
        Subtract of 2 maps
        '''
        if (type(other) is image2d):
            return image2d(self.field-other.field,self.res)
        elif (type(other) is float):
            return image2d(self.field-other,self.res)
    
    def __mul__(self,other):
        '''
        Multiply case by case
        '''
        if (type(other) is image2d):
            return image2d(self.field*other.field,self.res)
        if (type(other) is mask2d):
            return image2d(self.field*other.field,self.res)
        if (type(other) is float):
            return image2d(self.field*other,self.res)
        
    
    def __div__(self,other):
        'Divide self by other case by case'
        if (type(other) is image2d):
            return image2d(self.field*1/other.field,self.res)
        elif (type(other) in float):
            return self*1/other
    
    def pow(self, nb):
        '''
        map power nb
        
        :param nb:
        :type nb: float    
        '''
        
        return image2d(np.power(self.field,nb),self.res)
        
    def mask_build(self,polygone=False,r=0,grainId=[],pos_center=0):
        '''
        Create a mask map with NaN value where you don't want data and one were you want
        The default mask is a circle of center you choose and radius you define. 
        
        :param polygone: make a polygone mask ('not yet implemented')
        :type polygone: bool
        :param r: radius of the circle (warning what is the dimention of r mm ?)
        :type r: float
        :param grainId: You select the grainId you want in an array
        :type: array
        :return: mask
        :rtype: image2d
        :return: vec (vector of grainId is selction by grain or pos_center if selection by circle or 0 if polygone )
        :rtype: array
        
        .. note:: if you want applied a mask one your data just do data*mask where data is an image2d object
        '''
        # size of the figure
        ss=np.shape(self.field)
        mask_map=np.empty(ss, float)
        mask_map.fill(np.nan)
        
        # option 1 : draw polygone
        if polygone:
            print('not yet implemented')           
            xp=0
                    
        # option 2 : you want are circle
        elif r!=0:
            if np.size(pos_center)==1:
                self.plot()
                plt.waitforbuttonpress()
                print('clic to the center of the circle')
                xp=np.int32(np.array(plt.ginput(1))/self.res)
            else:
                xp=pos_center
            
            idx=[]
            plt.close('all')
            for i in np.int32(np.arange(2*r/self.res+1)+xp[0][0]-r/self.res):
                for j in np.int32(np.arange(2*r/self.res+1)+xp[0][1]-r/self.res):
                    if (((i-xp[0][0])**2+(j-xp[0][1])**2)**(0.5)<r/self.res):
                        idx.append([i,j])
            idx2=np.array(idx)
            y=ss[0]-idx2[:,1]
            x=idx2[:,0]
            v=(y>=0)*(y<ss[0])*(x>=0)*(x<ss[1])
            mask_map[[y[v],x[v]]]=1
            
            
            pc=float(sum(v))/float(len(y))
            if pc<1:
                print('WARNING : area is close to the border only '+str(pc)+'% of the initial area as been selected')
            
        # option 3 : grainId    
        else:
            if len(grainId)!=0:    
                gId=grainId
            else:
                plt.imshow(self.field,aspect='equal')
                plt.waitforbuttonpress()
                print('Select grains :')
                print('midle mouse clic when you are finish')
                xp=np.int32(np.array(plt.ginput(0)))
                plt.close('all')
                gId=self.field[xp[:,1],xp[:,0]]
            
            xp=gId
            for i in range(len(gId)):
                idx=np.where(self.field==gId[i])
                mask_map[idx]=1
        
        return mask2d(mask_map,self.res),xp
    
    
    def skeleton(self):
        '''
            Skeletonized a label map build by grain_label
        '''
        import micro2d
        # derived the image 
        a=self.diff('x')
        b=self.diff('y')
        # Normelized to one
        a=a/a
        b=b/b
        # Replace NaN by 0
        a.field[np.isnan(a.field)]=0
        b.field[np.isnan(b.field)]=0
        # Build the skeleton
        skel=a+b
        id=np.where(skel.field>0)
        skel.field[id]=1
        
        return micro2d(skel.field,skel.res)
        
    def vtk_export(self,nameId):
        '''
            Export the image2d into vtk file
            :param nameId: name of the output file
            :type name: str
        '''

        # size of the map
        ss=np.shape(self.field)
        # open micro.vtk file
        micro_out=open(nameId+'.vtk','w')
        # write the header of the file
        micro_out.write('# vtk DataFile Version 3.0 ' + str(datetime.date.today()) + '\n')
        micro_out.write('craft output \n')
        micro_out.write('ASCII \n')
        micro_out.write('DATASET STRUCTURED_POINTS \n')
        micro_out.write('DIMENSIONS ' + str(ss[1]) + ' ' + str(ss[0]) +  ' 1\n')
        micro_out.write('ORIGIN 0.000000 0.000000 0.000000 \n')
        micro_out.write('SPACING ' + str(self.res) + ' ' + str(self.res) + ' 1.000000 \n')
        micro_out.write('POINT_DATA ' + str(ss[0]*ss[1]) + '\n')
        micro_out.write('SCALARS scalars float \n')
        micro_out.write('LOOKUP_TABLE default \n')
        for i in list(range(ss[0]))[::-1]:
            for j in list(range(ss[1])):
                micro_out.write(str(int(self.field[i][j]))+' ')
            micro_out.write('\n')
        
                
        micro_out.close()
        
        return "vtk file created"
    
    def auto_correlation(self,pad=1):
        '''
        Compute the auto-correlation function from Dumoulin 2003
        
        :param pad: padding with mean data (default 1)
        :type pad: int
        
        :return: res_auto which the auto correlation function
        :rtype: image2d
        :return: Cinf see equation 6 Doumalin 2003
        :rtype: float
        :return: profil, it correspond to 180 profile taken for along line with direction between 0 and 180°. The profil[i,:] correspond to a line making an angle of i degree from x direction.
        :rtype: np.array 180*length profile
        :return: xi, it correspond to the distance along the profile for the variable profile
        :rtype: np.array 180*length profile
        :return: cross, it correspond to the length of the autocorrelation radius
        :rtype: np.array dim=180
        :Exemple:
            >>> # Using a fake image
            >>> data=image2d.image2d(np.eye(200),1)
            >>> [res_auto,Cinf,profil,xi,cross]=data.auto_correlation(pad=2)
            >>> # plot the auto correlation function
            >>> res_auto.plot()
            >>>plt.show()
            >>> # plot the auto_correlation radius in function of the angle of the profile
            >>> plt.plot(cross)
            >>> plt.xlabel('angle degree')
            >>> plt.ylabel('auto correlation radius (mm)')
            >>> plt.show()
            >>> # extract the max and min auto correlation radius
            >>> id=np.where(cross==np.min(cross))
            >>> minC=id[0][0]
            >>> id=np.where(cross==np.max(cross))
            >>> maxC=id[0][0]
            >>> ss=np.shape(xi) 
            >>> CC=np.ones(ss[1])*Cinf
            >>> plt.figure()
            >>> plt.plot(xi[0,:],CC)
            >>> plt.hold('on')
            >>> plt.plot(xi[minC,:],profil[minC,:])
            >>> plt.plot(xi[maxC,:],profil[maxC,:])
            >>> plt.show()
            
        '''
        # build a square picture nm*nm
        nm=np.min(np.shape(self.field))
        data=self.field[0:nm-1,0:nm-1]
        # mean of the data field
        mean_data=np.nanmean(data)
        # if there is nan value in data you need to replace them the choise now is to replace them by the mean value but may be interpolation can be better.
        id=np.isnan(data)
        if np.size(id)>0:
            data[id]=mean_data
        #
        fpad=np.ones([pad*nm,pad*nm])*mean_data
        fpad[0:nm-1,0:nm-1]=data
        ##
        FFT_fpad=np.fft.fft2(fpad)
        abs_FFTpad=np.abs(FFT_fpad)**2
        An=np.fft.ifft2(abs_FFTpad)
        mAn=np.nanmax(An)
        
        Autocor=np.abs(np.fft.fftshift(An/mAn));
        Cinf=mean_data**2/np.mean(fpad**2);
        res_auto=image2d(Autocor,1)
        
        
        # etract min an max direction
        ss=np.shape(res_auto.field)
        x0=ss[0]/2
        y0=x0
        cross=np.ones(180)
        for i in list(range(180)):
            xt=x0+np.cos(i*math.pi/180.)*x0
            yt=y0+np.sin(i*math.pi/180.)*y0
            pos=np.array([[x0,y0],[xt,yt]])
            [x,out,pos]=res_auto.extract_profil(pos=pos)
            if i==0:
                l=np.size(x)-1
                profil=np.ones([180,l])
                xi=np.ones([180,l])
            profil[i,:]=out[0:l]
            xi[i,:]=x[0:l]
            #
            id=np.where(profil[i,:]-Cinf<0)
            if np.size(id)==0:
                cross[i]=np.size(x)
            else:
                cross[i]=id[0][0]
            
        return res_auto,Cinf,profil,xi*self.res,cross*self.res
    def extract_profil(self,pos=0):
        '''
        Extract data along a given profile
        
        :param pos: correspond to the position of the line np.array([[X0,Y0],[X1,Y1]])
        :type pos: np.array dim 2*2
        :return: x
        :rtype: postion along the profile
        :return: out
        :rtype: value along the profile
        :return: pos
        :rtype: starting and ending points
        '''
        
        # size of the map
        ss=np.shape(self.field)
        # plot the data with phi1 value
        if np.size(pos)==1:
            h=plt.figure()
            self.plot()
            # select initial and final points for the line
            print('Select initial and final points for the line :')
            pos=np.array(pylab.ginput(2))
            plt.close(h)
        
        yy=np.float32([pos[0][0],pos[1][0]])/self.res
        xx=np.float32([pos[0][1],pos[1][1]])/self.res
        
        # numbers of pixel along the line
        nb_pixel=np.int32(np.sqrt((xx[1]-xx[0])**2+(yy[1]-yy[0])**2))
        
        # calcul for each pixel
        out=[]
        x=[]
        xi=[]
        yi=[]
        for i in list(range(nb_pixel)):
            # find the coordinate x an y along the line
            xi.append(ss[0]-np.int32(np.round(i*(xx[1]-xx[0])/nb_pixel+xx[0])))
            yi.append(np.int32(np.round(i*(yy[1]-yy[0])/nb_pixel+yy[0])))
            # extract phi and phi1
            out.append(self.field[xi[i],yi[i]])
            if i>0:
                x.append(np.sqrt((xi[i]-xi[0])**2+(yi[i]-yi[0])**2))
            else:
                x.append(0.0)
            
        return x, out, pos
        
        
        
###########################################
################# mask2d ##################
###########################################        
        
        
class mask2d(image2d):
    '''
    mask2d is a class which herit from image2d but it is restricted to microstructure object (background NaN, selected area 1)
    '''
    pass

    def __init__(self,field, resolution):
        '''
        Constructor : field is a matrix and resolution is the step size in mm
        
        :param field: tabular of scalar data
        :type field: array
        :param resolution: step size resolution (millimeters) 
        :type resolution: float
        '''
        
        if np.size(np.where((~np.isnan(field)) & (field!=1)))==0:
            self.field=field
            self.res=resolution
        else:
            print('error not only NaN and 1 in field')         
        return
    

    def load_mask(adr_bmp,res=1):
        '''
        Create a mask from a black and white bmp image where whith correspond to the selected area

        :param adr_bmp: path of the mask bmp file
        :type adr_bmp: str
        :param res: resolution of the picture mm (default 1)
        :type res: float
        :return mask:
        :rtype mask: mask2d 
        '''
        # Load bmp file
        image_bmp = io.imread(adr_bmp)

        tmp=image_bmp[:,:,0]
        # create nan matrix
        ss=np.shape(tmp)
        mask=np.zeros(ss)
        mask[:]=np.NaN
        # replace by one the white area
        id=np.where(tmp!=0)
        mask[id]=1

        return mask2d(mask,res)

    def complementary(self):
        '''
        Return complementary mask
        '''
        m=np.ones(np.shape(self.field))
        m[np.where(self.field==1)]=np.NaN

        return mask2d(m,self.res)

###########################################
################# micro2d ##################
###########################################    

class micro2d(image2d):
    '''
    micro2d is a class which herit from image2d but it is restricted to microstructure object (background 0, boundary 1)
    '''
    pass

    def __init__(self,field, resolution):
        '''
        Constructor : field is a matrix and resolution is the step size in mm
        
        :param field: tabular of scalar data
        :type field: array
        :param resolution: step size resolution (millimeters) 
        :type resolution: float
        '''
        
        if np.size(np.where((field!=0) & (field!=1)))==0:
            self.field=field
            self.res=resolution
        else:
            print('error not only 0 and 1 in field')         
        return
    
    
    def grain_label(self):
        '''
        Label area in a black and white picture
        
        .. note:: black color for the background and white color for the boundary
        '''
        # function which label a microstructure skeleton in one number per grain
        new_img=self.field
        res=self.res
        new_grain = morphology.label(new_img, connectivity=1, background=1)
        grains=image2d(new_grain,res)
        return grains
    
    def plotBoundary(self,dilatation=0):
        '''
        Add boundary to the figure
        
        :param dilatation: number of iteration for the dilation of 1 value - used this to make larger the boundaries (default 2)
        :type dilatation: int
        :Exemple:
            >>> data.phi1.plot()
            >>> data.micro.plotBoundary(dilatation=10)
            >>> plt.show()
        
        .. note:: plot only the value of the pixel equal at 1
        '''
        # extract microstructure matrix
        micro=self.field
        # make the dilation the number of time wanted
        if dilatation>0:
            micro=scipy.ndimage.binary_dilation(micro, iterations=dilatation)
        # create a mask with the 0 value
        micro = np.ma.masked_where(micro == 0, micro)
        # size of the image2d
        ss=np.shape(self.field)    
        # plot the boundaries
        plt.imshow(micro, extent=(0,ss[1]*self.res,0,ss[0]*self.res), interpolation='none',cmap=cm.gist_gray)
        
        return
    
    def area(self):
        '''
        Compute the grain area for each grain
        
        :return: g_arean array of scalar of grain area in mm^2
        :rtype: g_area np.array
        '''
        
        g_map=self.grain_label()
        
        g_area=np.zeros(np.max(g_map.field))
        
        for i in list(range(np.size(g_area))):
            g_area[i]=np.size(np.where(g_map.field==i))*g_map.res**2.
            
        return g_area

