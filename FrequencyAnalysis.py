'''
Created on 18 apr 2018

@author: jimmijamma
'''
import matplotlib as mpl 
mpl.use('TkAgg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import figure
import os
import numpy as np
from Wrapper import Wrapper
from scipy import interpolate
from datetime import datetime
import time
from ImgProc import ImgProc,twoD_Gaussian

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
        self.improc=ImgProc()
        
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
        
        self.display_fft_periods(periods,pos_spectrum,dir_name,param_name=y_name,peak_threshold=peak_threshold)
        self.display_fft_freqs(pos_freq,pos_spectrum,Ts,dir_name,param_name=y_name,peak_threshold=peak_threshold)
        
        return pos_freq,pos_spectrum


    def display_fft_periods(self,periods,spectrum,dir_name,param_name,peak_threshold=None):
        
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
        hours=(obs.timestamp-self.wrapper.observations[0].timestamp)*1.0/(60*60*1000*1000*1000)
        area_roi,aspect_roi,R,third_moment,uniformity=self.improc.ROI_analyseMoments(ROI)
        fp.write('%d %.6f %d %d %.6f %.6f %.6f %.6f\n' % (obs.id,hours,obs.ACrate,area_roi,aspect_roi,R,third_moment,uniformity))
        
    def writeROIcooccurrence(self,ROI,obs,fp):
        hours=(obs.timestamp-self.wrapper.observations[0].timestamp)*1.0/(60*60*1000*1000*1000)
        contrast, dissimilarity, homogeneity, energy, correlation, ASM=self.improc.ROI_analyseCooccurrenceMatrix(ROI,distances=[1],angles=[0])
        fp.write('%d %.6f %d %6f %.6f %.6f %.6f %.6f %.6f\n' % (obs.id,hours,obs.ACrate,contrast, dissimilarity, homogeneity, energy, correlation, ASM))
        
    def experiment(self,collection,x_dim,y_dim):
        txt_file_cooccurrence = open(self.txt_cooccurrence,'w')
        txt_file_moments = open(self.txt_moments,'w')
        
        for ii,o in enumerate(collection):   
            print "Resizing images, analysing and writing on disk: %d of %d" % (o.id,len(collection))  
            # original image
            image=o.window
            # interpolated image
            interp_image=self.improc.img_interpolateImage(image, x_dim, y_dim)
            # detecting blob
            ROI, row_interval, col_interval=self.improc.img_detectROI(interp_image, thr_factor=0.5, zero_padding=False)
            
            self.writeROImoments(ROI,o,txt_file_moments)
            self.writeROIcooccurrence(ROI,o,txt_file_cooccurrence)
            #print time.time()-t0
        
        txt_file_cooccurrence.close()   
        txt_file_moments.close()
        
        print "Experiment successfully completed!"
        print
            
          
    def experiment_with_resize_PCA(self, collection, x_dim, y_dim): 
        update_guess=10 # specify the n. of iteration after which start to update the guess, otherwise set to False
        guess=[20000,x_dim/2,y_dim/2,x_dim,y_dim,0,0]
        list_popt=[]
        if not os.path.exists(self.dir_name+'/PCAresults'):
            os.makedirs(self.dir_name+'/PCAresults')
        txt_file_cooccurrence = open(self.txt_cooccurrence,'w')
        txt_file_moments = open(self.txt_moments,'w')
        
        for ii, o in enumerate(collection):
            print "Resizing images, analysing and writing on disk: %d of %d" % (o.id,len(collection))
            t0 = time.time()
               
            img=o.window
            # interpolating image
            interp_img=self.improc.img_interpolateImage(img, x_dim, y_dim)
            # fitting the interpolated image wiht a 2D Gaussian function
            popt=self.improc.img_fitGaussian(img=interp_img, func=twoD_Gaussian, guess=guess)
            # popt = (amplitude,dx,dy,sigma_x,sigma_y,theta,offset)
            list_popt.append(popt)
            
            mu_x=popt[1]
            mu_y=popt[2]
            s_x=popt[3]**2
            s_y=popt[4]**2
    
            # resizing the profile using PCA
            resized_img=self.improc.img_resizeProfilePCA(interp_img, mu_x, mu_y, s_x, s_y, plot_results=True, it=ii, folder=self.dir_name+'/PCAresults')
     
            if update_guess!=False:
                if ii>update_guess:
                    guess=np.mean(list_popt,axis=0)
                
            res=self.improc.img_detectROI(resized_img, thr_factor=0.5, zero_padding=False)
            if res==-1:
                print "ERROR detecting ROI!!"
                continue
            
            ROI,row_interval,col_interval=res
            
            self.writeROImoments(ROI,o,txt_file_moments)
            self.writeROIcooccurrence(ROI,o,txt_file_cooccurrence)
            #print time.time()-t0
        
        txt_file_cooccurrence.close()   
        txt_file_moments.close()
    
        print "Experiment successfully completed!"
        print  
        
        
    def experiment_with_resize_aspect(self, collection, x_dim, y_dim):
        
        txt_file_cooccurrence = open(self.txt_cooccurrence,'w')
        txt_file_moments = open(self.txt_moments,'w')
        
        for ii,o in enumerate(collection):
            print "Resizing images, analysing and writing on disk: %d of %d" % (ii,len(collection))
            interp_image=self.improc.img_interpolateImage(o.window,x_dim,y_dim)
            res=self.improc.img_detectROI(interp_image, thr_factor=0.5, zero_padding=False)
            if res==-1:
                print "ERROR!!"
                continue
            ROI,row_interval,col_interval=res
            roi=np.array(ROI)
            height=np.shape(roi)[0]
            width=np.shape(roi)[1]
            aspect_roi=1.0*width/height
            # resizing image
            resized_img=self.improc.img_resizeProfileROIaspect(o.window, aspect_roi, x_dim=x_dim, y_dim_=y_dim)
            # detecting image ROI
            res=self.improc.img_detectROI(resized_img, thr_factor=0.5, zero_padding=False)
            if res==-1:
                print "ERROR!!"
                continue
            ROI,row_interval,col_interval=res
         
            self.writeROImoments(ROI, o, txt_file_moments)
            self.writeROIcooccurrence(ROI, o, txt_file_cooccurrence)
         
        txt_file_cooccurrence.close()   
        txt_file_moments.close()
        
        print "Experiment successfully completed!"
        print 
        
        