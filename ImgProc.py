'''
Created on 18 mag 2018

@author: jimmijamma
'''
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
import cv2
from skimage.feature import greycomatrix, greycoprops

class ImgProc(object):
    '''
    Class that implements method for processing single images
    '''

    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    
    def img_interpolateImage(self,image,x_dim=1800,y_dim=1200, algorithm=cv2.INTER_LANCZOS4):
        '''
        Method for resizing the image using interpolation
        '''
        interp_image=cv2.resize(image,(x_dim,y_dim),interpolation=algorithm)
        return interp_image
    
    
    def img_detectROI(self,image,thr_factor=0.5,zero_padding=False):
        '''
        Method for detecting the PSF profiles defining a ROI
        '''
        im=np.array(image)
        fwhm_list=[] # points of the ROI. For each row: (row_coordinate,first_column,last_column)
        idy=0
        threshold=np.nanmax(im)*thr_factor # FWHM threshold
        for row in im:
            idx=0
            flagx=False
            flagy=False
            first=second=None
            for pix in row:           
                if pix>threshold and flagx==False:
                    first=idx
                    flagx=True
                if pix<threshold and flagx==True:
                    second=idx-1
                    flagx=False
                    flagy=True
                    break
                idx+=1                
            if first!=None and second!=None:
                fwhm_list.append([idy,first,second])
            idy+=1
            
        # defining ROI
        rows_list=[row[0] for row in fwhm_list]
        row_interval=[np.min(rows_list),np.max(rows_list)]
        col_interval=[np.min([row[1] for row in fwhm_list]),np.max([row[2] for row in fwhm_list])]
        
        ROI=im[row_interval[0]:row_interval[1],col_interval[0]:col_interval[1]]
        if zero_padding==True:
            low_values_flags = ROI < threshold  # Where values are low
            ROI[low_values_flags] = 0
          
        return ROI, row_interval, col_interval
    
    
    def img_resizeProfileROIaspect(self,image,aspectROI,algorithm=cv2.INTER_LANCZOS4):
        '''
        Method for shrinking PSF profiles through ROI aspect
        '''
        x_dim=1800
        y_dim=int(1200*aspectROI)
        interp_image=cv2.resize(image,(x_dim,y_dim),interpolation=algorithm)
        return interp_image
    
    
    def img_resizeProfilePCA(self,image,mu_x,mu_y,s_x,s_y,plot_results=False,it=1,folder=''):
        '''
        Method for shrinking PSF profiles using PCA
        '''
        indices=np.where(image)
        X=indices[1]
        Y=indices[0]
        
        #corr=np.corrcoef(X, Y)
        #corr=corr[1,0]
        #cov= [[s_x,corr*np.sqrt(s_x)*np.sqrt(s_y)],[corr*np.sqrt(s_x)*np.sqrt(s_y),s_y]]
        cov= [[s_x,0],[0,s_y]] # covariance matrix
        L,U=np.linalg.eigh(cov)
        z1=U[0,0]*(X-mu_x)+U[0,1]*(Y-mu_y) 
        z2=U[1,0]*(X-mu_x)+U[1,1]*(Y-mu_y)
        
        L[1]=L[1]/L[0]
        L[0]=1
        
        w1,w2=(z1/np.sqrt(L[0]))+mu_x,(z2/np.sqrt(L[1]))+mu_y # new sets of coordinates 
        
        x_lim=np.shape(image)[1]
        y_lim=np.shape(image)[0]
        x = np.linspace(start=0, stop=x_lim-1, num=x_lim)
        y = np.linspace(start=0, stop=y_lim-1, num=y_lim)
        x, y = np.meshgrid(x, y)
        
        resized_img=interpolate.griddata((w1,w2), image.ravel(), xi=(x,y), method='cubic', fill_value=np.min(image))
        
        if plot_results==True:
            fig, axes = plt.subplots(nrows=1, ncols=2)
            axes[0].imshow(image)
            axes[0].set_title('Original (interpolated)')
            axes[1].imshow(resized_img)
            axes[1].set_title('Shrinked with PCA')
            fig.suptitle('Image #%i' %it)
            plt.savefig(folder+'/PCAresult'+str(it)+'.png')
            plt.close()
        
        return resized_img
    
        
    def ROI_analyseMoments(self,ROI):
        '''
        Method for analysing moments and shape features of the ROI
        '''
        roi=np.array(ROI).astype(int).astype(float)
        height=np.shape(roi)[0]
        width=np.shape(roi)[1]
        area_roi=height*width
        aspect_roi=1.0*width/height
        L=65536 # number of levels for 16-bit
        hist,bins = np.histogram(roi.ravel(),L,[0,L])
        norm_hist=[]
        hist[0]=0 # non considero i pixel a 0
        s=sum(hist)
        for el in hist:
            if el>0:
                norm_hist.append(el*1.0/s)
            else:
                norm_hist.append(0.0)
                
        roi[roi==0]=np.nan
        m=np.nanmean(roi.ravel())
        
        second_moment=0
        third_moment=0
        uniformity=0
        entropy=0
        for i in range(L):
            second_moment+=((i-m)**2*norm_hist[i])
            third_moment+=((i-m)**3*norm_hist[i])
            uniformity+=norm_hist[i]**2
            #entropy-=norm_hist[i]*np.log2(norm_hist[i])
            
        second_moment=second_moment/((L-1)**2)   
        R=1-1.0/(1+second_moment) # relative smoothness
        
        return area_roi,aspect_roi,R,third_moment,uniformity
    
    
    def ROI_analyseCooccurrenceMatrix(self,ROI,distances,angles):
        '''
        Method for analysing co-occurrence matrix properties of the ROI.
        In this case the ROI is evaluated as a 8-bit image
        '''
        glcm = greycomatrix(ROI.astype(int)/256, distances, angles, 256, symmetric=True, normed=True)
        contrast = greycoprops(glcm, prop='contrast')[0, 0]
        dissimilarity = greycoprops(glcm, prop='dissimilarity')[0, 0]
        homogeneity = greycoprops(glcm, prop='homogeneity')[0, 0]
        energy = greycoprops(glcm, prop='energy')[0, 0]
        correlation = greycoprops(glcm, prop='correlation')[0, 0]
        ASM = greycoprops(glcm, prop='ASM')[0, 0] 
        return contrast,dissimilarity,homogeneity,energy,correlation,ASM
    
    
        