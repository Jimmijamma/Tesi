'''
Created on 18 apr 2018

@author: jimmijamma
'''
import matplotlib as mpl 
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import figure
import os
import numpy as np
from skimage.feature import greycomatrix, greycoprops
import cv2
from Wrapper import Wrapper
import json
from scipy import interpolate

class FrequencyAnalysis(object):
    '''
    classdocs
    '''
    def __init__(self, wrapper):
        '''
        Constructor
        '''
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
    
    def writeResultsMoments(self,n_elements=None):
        collection=self.wrapper.getCollection(n_elements)
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
        print
        
    def writeResultsCooccurrence(self,n_elements=None):
        collection=self.wrapper.getCollection(n_elements)
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
        print
        
    def interpolate_measurements(self,xdata,ydata):
        new_xdata = np.linspace(min(xdata), max(xdata), len(xdata))
        f = interpolate.interp1d(xdata, ydata)
        new_ydata = f(new_xdata)   # use interpolation function returned by `interp1d`
        return new_xdata,new_ydata
    
    def evaluate_fft(self,x_data,y_data,y_name,peak_threshold,dir_name):
    
        Ts=(max(x_data)-min(x_data))/len(x_data)
    
        A=np.fft.fft(y_data)
        freq = np.fft.fftfreq(len(x_data))
        spectrum=np.abs(A)
        
        pos_freq=freq[:len(freq)/2]
        pos_spectrum=spectrum[:len(spectrum)/2]
        
        periods=[]
        for f in pos_freq:
                periods.append((Ts/f))
        
        threshold = peak_threshold * max(pos_spectrum)
        
        #mask = pos_spectrum > threshold 
        #peaks = periods[mask]
        
        self.display_fft(periods,pos_spectrum,dir_name,param_name=y_name,peak_threshold=peak_threshold)

        return pos_freq,pos_spectrum


    def display_fft(self,periods,spectrum,dir_name,param_name,peak_threshold=None):
        
        fig,ax=plt.subplots()
        ax.plot(periods,spectrum, color='green')
        ax.set_xlabel('Period [h]')
        ax.set_ylabel('Amplitude')
        plt.ylim([None, max(spectrum)*1.1])
        
        bbox_props = dict(boxstyle='round', alpha=0.7, ec="b", lw=1, fc='white')
        max_x=np.max(periods)
        max_y=np.max(spectrum)
        top_peaks=np.sort(spectrum)[-3:]
        for i,s in enumerate(spectrum):
            if s>=0.05*max(spectrum) and s in top_peaks:
                ax.text(periods[i], spectrum[i], "T="+str(round(periods[i],3)), ha="center", va="center",
                size=10, bbox=bbox_props)
        
        plt.title("Periodicity of '%s' parameter" % param_name)
        plt.savefig(dir_name+'/fft_'+param_name+'.png')
        plt.close()
    
    def display_timedomain(self,x_var,y_var,y_name,ACmotion, dir_name):
        
        w, h = figure.figaspect(6.)
        fig, ax = plt.subplots(figsize=(h,w))
    
        # Twin the x-axis twice to make independent y-axes.
        axes = [ax, ax.twinx()]
        
        # Make some space on the right side for the extra y-axis.
        fig.subplots_adjust(right=0.75)
        axes[0].set_xlabel('hours')
        axes[0].plot(x_var,y_var, color='green')
        axes[0].tick_params(axis='y', colors='green')
        axes[0].set_ylabel(y_name, color='green')
        
        axes[1].plot(x_var,ACmotion, color='red')
        axes[1].tick_params(axis='y', colors='red')
        axes[1].set_ylabel('AC motion',color='red')
        
        plt.savefig(dir_name+'/ACmotionVS'+y_name+'.png')
        plt.close()
        
    def readResultsMoments(self):
        print "Reading and displaying moments results..."
        f=open(self.dir_name+'/resultsMoments.txt','r')
        measure_list=[]
        for line in f:
            l=line.split()
            m=[l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]]
            measure_list.append(m)
        f.close()
        
        dir_name=self.dir_name+'/moments'
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        
        ids=[int(row[0]) for row in measure_list]
        hours=[float(row[1]) for row in measure_list]
        ACmotion=[int(row[2]) for row in measure_list]
        area_roi=[int(row[3]) for row in measure_list]
        aspect_roi=[float(row[4]) for row in measure_list]
        R=[float(row[5]) for row in measure_list]
        third_moment=[float(row[6]) for row in measure_list]
        uniformity=[float(row[7]) for row in measure_list]
        
        vars_dict=dict({'area_roi':area_roi,'aspect_roi':aspect_roi,'R':R,'third_moment':third_moment,'uniformity':uniformity,})
        for v in vars_dict:
            xvar=hours
            yvar=vars_dict[v]
            interp_x,interp_y=self.interpolate_measurements(xvar, yvar-np.mean(yvar))
            self.display_timedomain(interp_x,interp_y,v,ACmotion,dir_name)
            freq,spectrum=self.evaluate_fft(interp_x,interp_y,v,0.05,dir_name)
            
        print "Results read and displayed successfully!"
        print
            
    def readResultsCooccurrence(self):
        print "Reading and displaying co-occurrence results..."
        f=open(self.dir_name+'/resultsCooccurrenceMatrix.txt','r')
        measure_list=[]
        for line in f:
            l=line.split()
            m=[l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8]]
            measure_list.append(m)
        f.close()
        
        dir_name=self.dir_name+'/cooccurrence_matrix'
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        
        ids=[int(row[0]) for row in measure_list]
        hours=[float(row[1]) for row in measure_list]
        ACmotion=[int(row[2]) for row in measure_list]
        contrast=[float(row[3]) for row in measure_list]
        dissimilarity=[float(row[4]) for row in measure_list]
        homogeneity=[float(row[5]) for row in measure_list]
        energy=[float(row[6]) for row in measure_list]
        correlation=[float(row[7]) for row in measure_list]
        ASM=[float(row[7]) for row in measure_list]
        
        vars_dict=dict({'contrast':contrast,'dissimilarity':dissimilarity,'homogeneity':homogeneity,'energy':energy,'correlation':correlation,'ASM':ASM})
        for v in vars_dict:
            xvar=hours
            yvar=vars_dict[v]
            interp_x,interp_y=self.interpolate_measurements(xvar, yvar-np.mean(yvar))
            self.display_timedomain(interp_x,interp_y,v,ACmotion,dir_name)
            freq,spectrum=self.evaluate_fft(interp_x,interp_y,v,0.05,dir_name)
        
        print "Results read and displayed successfully!"
        print
    
if __name__ == '__main__':
    json_data=open('CalWrapper.json')
    json_obj = json.load(json_data)
    json_data.close()
    
    print
    w=Wrapper(json_obj)
    freq_analysis=FrequencyAnalysis(w)
    freq_analysis.writeResultsMoments(10)
    freq_analysis.writeResultsCooccurrence(10)
    freq_analysis.readResultsCooccurrence()
    freq_analysis.readResultsMoments()
    
    