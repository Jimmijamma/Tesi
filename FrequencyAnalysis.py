'''
Created on 18 apr 2018

@author: jimmijamma
'''

import numpy as np
from skimage.feature import greycomatrix, greycoprops
import cv2
from Wrapper import Wrapper
import json

class FrequencyAnalysis(object):
    '''
    classdocs
    '''
    def __init__(self, wrapper):
        '''
        Constructor
        '''
        import os
        self.wrapper=wrapper
        self.dir_name=self.wrapper.dir_name+'/FrequencyAnalysis'
        if not os.path.exists(self.dir_name):
            os.makedirs(self.dir_name)
        
    def interpolateImage(self,image,x_dim=1800,y_dim=1200):
        interp_image=cv2.resize(image,(x_dim,y_dim),interpolation=cv2.INTER_LANCZOS4)
        return interp_image
        
    def detectROI(self,image,zero_padding=True):
        im=np.array(image)
        fwhm_list=[] # points of the ROI. For each row: (row_coordinate,first_column,last_column)
        idy=0
        threshold=np.max(im)*1.0/2 # FWHM threshold
        for row in im:
            idx=0
            flagx=False
            flagy=False
            first=second=None
            for pix in row:           
                if pix>=threshold and flagx==False:
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
        row_interval=[np.min([row[0] for row in fwhm_list]),np.max([row[0] for row in fwhm_list])]
        col_interval=[np.min([row[1] for row in fwhm_list]),np.max([row[2] for row in fwhm_list])]
        
        ROI=im[row_interval[0]:row_interval[1],col_interval[0]:col_interval[1]]
        if zero_padding==True:
            low_values_flags = ROI < threshold  # Where values are low
            ROI[low_values_flags] = 0
          
        return ROI, row_interval, col_interval
    
    def analyseROImoments(self,ROI):
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
    
    def analyseROIcooccurrence(self,ROI,distances,angles):
        glcm = greycomatrix(ROI.astype(int)/256, distances, angles, 256, symmetric=True, normed=True)
        contrast = greycoprops(glcm, prop='contrast')[0, 0]
        dissimilarity = greycoprops(glcm, prop='dissimilarity')[0, 0]
        homogeneity = greycoprops(glcm, prop='homogeneity')[0, 0]
        energy = greycoprops(glcm, prop='energy')[0, 0]
        correlation = greycoprops(glcm, prop='correlation')[0, 0]
        ASM = greycoprops(glcm, prop='ASM')[0, 0]
        
        return contrast,dissimilarity,homogeneity,energy,correlation,ASM
    
    def writeResultsMoments(self):
        collection=self.wrapper.observations[:10]
        txt_file = open(self.dir_name+'/resultsMoments.txt','w') 
        for o in collection:     
            # original image
            image=o.window
            # interpolated image
            interp_image=self.interpolateImage(image)
            # detecting blob
            ROI, row_interval, col_interval=self.detectROI(interp_image)
            hours=(o.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
            # analysing ROI moments
            area_roi,aspect_roi,R,third_moment,uniformity=self.analyseROImoments(ROI)
            txt_file.write('%d %.6f %d %d %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACmotion,area_roi,aspect_roi,R,third_moment,uniformity))
    
            print "Writing results of image moments: %d of %d" % (o.id,len(collection))
        txt_file.close()
        print "Results recorded successfully!"
        
    def writeResultsCooccurrence(self):
        collection=self.wrapper.observations
        txt_file = open(self.dir_name+'/resultsCooccurrenceMatrix.txt','w') 
        for o in collection:     
            # original image
            image=o.window
            # interpolated image
            interp_image=self.interpolateImage(image)
            # detecting blob
            ROI, row_interval, col_interval=self.detectROI(interp_image)
            hours=(o.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
            # analysing ROI moments
            contrast, dissimilarity, homogeneity, energy, correlation, ASM=self.analyseROIcooccurrence(ROI,distances=[1],angles=[0])
            txt_file.write('%d %.6f %d %6f %.6f %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACmotion,contrast, dissimilarity, homogeneity, energy, correlation, ASM))
    
            print "Writing results of co-occurrence matrices: %d of %d" % (o.id,len(collection))
        txt_file.close()
        print "Results recorded successfully!"
        
    
if __name__ == '__main__':
    json_data=open('CalWrapper.json')
    json_obj = json.load(json_data)
    json_data.close()
    
    print
    w=Wrapper(json_obj)
    freq_analysis=FrequencyAnalysis(w)
    freq_analysis.writeResultsMoments()
    
    