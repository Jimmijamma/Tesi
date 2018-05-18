'''
Created on 18 apr 2018

@author: jimmijamma
'''
import matplotlib as mpl 
from numpy import interp
from _ctypes import resize
mpl.use('TkAgg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import figure
import os
import numpy as np
from skimage.feature import greycomatrix, greycoprops
import cv2
from Wrapper import Wrapper
import json
from scipy import interpolate
from datetime import datetime
from scipy import optimize as opt
import time

class FrequencyAnalysis(object):
    '''
    classdocs
    '''
    def __init__(self, wrapper):
        '''
        Constructor
        '''
        self.wrapper=wrapper
        self.date=self.computeDate()
        
        self.deleteEmptyFolders()   
        
        self.dir_name=self.wrapper.dir_name+'/FrequencyAnalysis/'+self.date
        if not os.path.exists(self.dir_name):
            os.makedirs(self.dir_name)
            
        self.txt_cooccurrence=self.dir_name+'/resultsCooccurrenceMatrix.txt'
        self.txt_moments=self.dir_name+'/resultsMoments.txt'
        
        
    def deleteEmptyFolders(self):
        if os.path.exists(self.wrapper.dir_name+'/FrequencyAnalysis/'):
            for x in os.walk(self.wrapper.dir_name+'/FrequencyAnalysis/'):
                if os.listdir(x[0])==[]:
                    if os.path.isdir(x[0]):
                        os.rmdir(x[0])
                        
    def computeDate(self):
        date_time=str(datetime.now()).split(' ')
        date=date_time[0]
        hour=date_time[1].split(':')[:-1]
        str_date=date+'_'+hour[0]+':'+hour[1]
        return str_date
        
    def interpolateImage(self,image,x_dim=1800,y_dim=1200, algorithm=cv2.INTER_LANCZOS4):
        interp_image=cv2.resize(image,(x_dim,y_dim),interpolation=algorithm)
        return interp_image
    
    def resizeProfile(self,image,aspectROI, algorithm=cv2.INTER_LANCZOS4):
        x_dim=1800
        y_dim=int(1200*aspectROI)
        interp_image=cv2.resize(image,(x_dim,y_dim),interpolation=algorithm)
        return interp_image
        
    def detectROI(self,image,thr_factor=0.5,zero_padding=False):
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
        
        pos_freq=freq[:len(freq)/2]/Ts
        pos_spectrum=spectrum[:len(spectrum)/2]
        
        periods=[]
        for f in pos_freq:
                periods.append((1.0/f))
        
        threshold = peak_threshold * max(pos_spectrum)
        
        #mask = pos_spectrum > threshold 
        #peaks = periods[mask]
        
        self.display_fft(periods,pos_spectrum,dir_name,param_name=y_name,peak_threshold=peak_threshold)
        self.display_fft_freqs(pos_freq,pos_spectrum,Ts,dir_name,param_name=y_name,peak_threshold=peak_threshold)
        
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
        
    def display_fft_freqs(self,freqs,spectrum,Ts,dir_name,param_name,peak_threshold=None):
        
        fig,ax=plt.subplots()
        ax.plot(freqs,spectrum, color='green')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Amplitude')
        plt.ylim([None, max(spectrum)*1.1])
        
        bbox_props = dict(boxstyle='round', alpha=0.7, ec="b", lw=1, fc='white')
        max_x=np.max(freqs)
        max_y=np.max(spectrum)
        top_peaks=np.sort(spectrum)[-3:]
        for i,s in enumerate(spectrum):
            if s>=0.05*max(spectrum) and s in top_peaks:
                ax.text(freqs[i], spectrum[i], "T="+str(round(1.0/(freqs[i]),3)), ha="center", va="center",
                size=10, bbox=bbox_props)
        
        plt.title("Periodicity of '%s' parameter" % param_name)
        plt.savefig(dir_name+'/fft_freqs_'+param_name+'.png')
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
        
    def readResultsMoments(self, timestamp=None):
        print "Reading and displaying moments results..."
        if timestamp!=None:
            f=open(self.wrapper.dir_name+'/FrequencyAnalysis/'+timestamp+'/resultsMoments.txt','r')
        else:
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
            self.display_timedomain(xvar,yvar,v,ACmotion,dir_name)
            freq,spectrum=self.evaluate_fft(interp_x,interp_y,v,0.05,dir_name)
            
        print "Results read and displayed successfully!"
        print
            
    def readResultsCooccurrence(self, timestamp=None):
        print "Reading and displaying co-occurrence results..."
        if timestamp!=None:
            f=open(self.wrapper.dir_name+'/FrequencyAnalysis/'+timestamp+'/resultsCooccurrenceMatrix.txt','r')
        else:
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
            self.display_timedomain(xvar,yvar,v,ACmotion,dir_name)
            freq,spectrum=self.evaluate_fft(interp_x,interp_y,v,0.05,dir_name)
        
        print "Results read and displayed successfully!"
        print
    
    
    def writeROImoments(self,ROI,obs,fp):
    #analysing ROI moments
        hours=(obs.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
        area_roi,aspect_roi,R,third_moment,uniformity=self.analyseROImoments(ROI)
        fp.write('%d %.6f %d %d %.6f %.6f %.6f %.6f\n' % (obs.id,hours,obs.ACrate,area_roi,aspect_roi,R,third_moment,uniformity))
        
    def writeROIcooccurrence(self,ROI,obs,fp):
        hours=(obs.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
        contrast, dissimilarity, homogeneity, energy, correlation, ASM=self.analyseROIcooccurrence(ROI,distances=[1],angles=[0])
        fp.write('%d %.6f %d %6f %.6f %.6f %.6f %.6f %.6f\n' % (obs.id,hours,obs.ACrate,contrast, dissimilarity, homogeneity, energy, correlation, ASM))
          
    def experiment_with_resize_PCA(self, collection, x_dim, y_dim): 
        update_guess=10 # specify the n. of iteration after which start to update the guess, otherwise set to False
        guess=[20000,x_dim/2,y_dim/2,x_dim,y_dim,0,0]
        list_popt=[]
        if not os.path.exists(self.dir_name+'/PCAresults'):
            os.makedirs(self.dir_name+'/PCAresults')
        txt_file_cooccurrence = open(self.txt_cooccurrence,'w')
        txt_file_moments = open(self.txt_moments,'w')
        
        for ii, o in enumerate(collection):
            t0 = time.time()
               
            img=o.window
            # interpolating image
            interp_img=self.interpolateImage(img, x_dim, y_dim)
            # fitting the interpolated image wiht a 2D Gaussian function
            popt=fitGaussian(img=interp_img, func=twoD_Gaussian, guess=guess)
            # popt = (amplitude,dx,dy,sigma_x,sigma_y,theta,offset)
            list_popt.append(popt)
            
            mu_x=popt[1]
            mu_y=popt[2]
            s_x=popt[3]**2
            s_y=popt[4]**2
    
            # resizing the profile using PCA
            resized_img=resizePCA(interp_img, mu_x, mu_y, s_x, s_y, plot_results=True, it=ii, folder=self.dir_name+'/PCAresults')
     
            if update_guess!=False:
                if ii>update_guess:
                    guess=np.mean(list_popt,axis=0)
                
            res=self.detectROI(resized_img, thr_factor=0.5, zero_padding=False)
            if res==-1:
                print "ERROR detecting ROI!!"
                continue
            
            ROI,row_interval,col_interval=res
            
            self.writeROImoments(ROI,o,txt_file_moments)
            self.writeROIcooccurrence(ROI,o,txt_file_cooccurrence)
            #print time.time()-t0
            print '%i of %i' %(ii+1,len(collection))
        
        txt_file_cooccurrence.close()   
        txt_file_moments.close()
    
        print "Experiment successfully completed!"
        print    
    
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def fitGaussian(img, func, guess):
    x_lim=np.shape(img)[1]
    y_lim=np.shape(img)[0]
    x = np.linspace(start=0, stop=x_lim-1, num=x_lim)
    y = np.linspace(start=0, stop=y_lim-1, num=y_lim)
    x, y = np.meshgrid(x, y)
    # constraints for fitting ((lower bounds),(upper bounds))
    bounds=((0,0,0,0,0,0,-np.inf),(80000,x_lim-1,y_lim-1,np.inf,np.inf,0.5,np.inf))
    popt, pcov = opt.curve_fit(func, (x, y), img.ravel(), p0=guess, bounds=bounds,maxfev=1000000)
    return popt

def resizePCA(img, mu_x, mu_y, s_x, s_y, plot_results=False, it=None, folder=None):
    indices=np.where(img)
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
    
    x_lim=np.shape(img)[1]
    y_lim=np.shape(img)[0]
    x = np.linspace(start=0, stop=x_lim-1, num=x_lim)
    y = np.linspace(start=0, stop=y_lim-1, num=y_lim)
    x, y = np.meshgrid(x, y)
    
    resized_img=interpolate.griddata((w1,w2), img.ravel(), xi=(x,y), method='cubic', fill_value=np.min(img))
    
    if plot_results==True:
        fig, axes = plt.subplots(nrows=1, ncols=2)
        axes[0].imshow(img)
        axes[0].set_title('Original (interpolated)')
        axes[1].imshow(resized_img)
        axes[1].set_title('Shrinked with PCA')
        fig.suptitle('Image #%i' %it)
        plt.savefig(folder+'/PCAresult'+str(it)+'.png')
        plt.close()
    
    return resized_img
    
if __name__ == '__main__':
    
    print
    w=Wrapper('CalWrapper.json')
    
    fa=FrequencyAnalysis(w)
    
    collection=w.getCollection()
    
    fa.experiment_with_resize_PCA(collection,x_dim=180,y_dim=120)
    
    '''
    for o in collection:
        print "Resizing images, analysing and writing on disk: %d of %d" % (o.id,len(collection))
        interp_image=freq_analysis.interpolateImage(o.window, 1800, 1200)
        res=freq_analysis.detectROI(interp_image, thr_factor=0.5, zero_padding=False)
        if res==-1:
            print "ERROR!!"
            continue
        ROI,row_interval,col_interval=res
        roi=np.array(ROI)
        height=np.shape(roi)[0]
        width=np.shape(roi)[1]
        aspect_roi=1.0*width/height
        # resizing image
        resized_img=freq_analysis.resizeProfile(o.window, aspect_roi)
        # detecting image ROI
        res=freq_analysis.detectROI(resized_img, thr_factor=0.5, zero_padding=False)
        if res==-1:
            print "ERROR!!"
            continue
        ROI,row_interval,col_interval=res
     
        hours=(o.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
        
        # analysin ROI moments
        txt_file_moments = open(freq_analysis.dir_name+'/resultsMoments.txt','a') 
        area_roi,aspect_roi,R,third_moment,uniformity=freq_analysis.analyseROImoments(ROI)
        if aspect_roi!=1.0:
            print "ERROR!"
            continue
        txt_file_moments.write('%d %.6f %d %d %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACrate,area_roi,aspect_roi,R,third_moment,uniformity))
        
        # analysing cooccurrence matrix
        txt_file_cooccurrence = open(freq_analysis.dir_name+'/resultsCooccurrenceMatrix.txt','a')
        contrast, dissimilarity, homogeneity, energy, correlation, ASM=freq_analysis.analyseROIcooccurrence(ROI,distances=[1],angles=[0])
        txt_file_cooccurrence.write('%d %.6f %d %6f %.6f %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACrate,contrast, dissimilarity, homogeneity, energy, correlation, ASM))
     
    txt_file_cooccurrence.close()   
    txt_file_moments.close()
    
    
    print "Finished resizing images!"
    print

    freq_analysis.readResultsCooccurrence()
    freq_analysis.readResultsMoments()
    
    '''

    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    
    

    
    
    
    