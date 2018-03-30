'''
Created on 23 gen 2018

@author: jimmijamma
'''
import numpy as np
import matplotlib as mpl 
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime

from Observation import Observation

pix_size_AL=47260/4500
pix_size_AC=60000/1966 # micron

class Wrapper(object):
    '''
    classdocs
    '''

    def __init__(self, json_obj):
        '''
        Constructor
        '''
        self.id = None
        self.first_obs = json_obj['first_transit'] # ns
        self.last_obs = json_obj['last_transit'] # ns
        self.time_interval = 1.0*(self.last_obs-self.first_obs)/(60*60*1000*1000*1000)
        self.n_transits = json_obj['n_transits']
        self.ccdRow = json_obj['ccdRow']
        self.ccdStrip = json_obj['ccdStrip']
        self.mag = json_obj['magnitude']
        self.wavn = json_obj['wavenumber']
        fov = json_obj['fov']
        if fov == True:
            self.fov=1
        else:
            self.fov=0
        self.tdi_period = None
        self.observations=self.parseJSON(json_obj)
        
        print "Wrapper loaded successfully"
        print "> CCD coordinates: (" + str(self.ccdRow) + "," + str(self.ccdStrip) + ")"
        print "> Field Of View (FOV): " +  str(self.fov)
        print "> Stars of magnitude interval " + str(self.mag) + " mag, wavenumber " + str(self.wavn) + " nm"
        print "> N. of transits: " + str(self.n_transits)
        print "> First transit on: " + str(self.first_obs)
        #datetime.datetime.fromtimestamp(self.first_obs/1000.0).strftime('%Y-%m-%d %H:%M:%S.%f')
        print "> Last transit on: " + str(self.last_obs)
        print "> Time interval: " + str(self.time_interval) + " hours"
        #datetime.datetime.fromtimestamp(self.last_obs/1000.0).strftime('%Y-%m-%d %H:%M:%S.%f')
        print
        
        
    def parseJSON(self, json_obj):
        obs_list=[]
        js_obs=json_obj['observations']
        for obs in js_obs:
            win=np.array(obs['window'])
            calCentroid_AC=obs['calCentroidAC']/pix_size_AC+5.5
            calCentroid_AL=obs['calCentroidAL']/pix_size_AL+8.5
            o = Observation(self, obs['id'], np.swapaxes(win,0,1), obs['gating'], obs['transitid'], obs['timestamp'], calCentroid_AC, calCentroid_AL)
            obs_list.append(o)
        return obs_list
        
    def createScatterPlot(self, path, filename):
        figg = plt.figure()
        axxx = figg.add_subplot(111, projection='3d')
        #axxx.scatter(dx.ravel(), dy.ravel(), values.ravel(),c=values.ravel(),cmap=plt.get_cmap('autumn'))
        figg.savefig(filename)
        
    def create2Dhistogram(self, collection, filename, dx,dy,values,xbins, ybins, normed):
        fig,ax=plt.subplots()
        bins=[np.linspace(-10,9,xbins),np.linspace(-7, 7,ybins)]
        my_cmap = plt.get_cmap('YlOrRd')
        my_cmap.set_under('w')
        h, xe, ye, im= ax.hist2d(dx.ravel(), dy.ravel(), bins=bins, cmap=my_cmap, weights=values.ravel(), cmin=0.1)
        ax.axhline(xmin=0.95,color='black')
        ax.axhline(xmax=0.05,color='black')
        ax.axvline(ymin=0.95,color='black')
        ax.axvline(ymax=0.05,color='black')
        plt.colorbar(im, ax=ax)
        plt.xlabel('AL direction')
        plt.ylabel('AC direction')
        plt.title('Stacked version of %s transits, with resolution %sx%s'%(len(collection),ybins,xbins))
        fig.savefig(filename)
        return xe, ye, h
    
    def estimateCentroids_fit(self, collection):
        c_list_x=[]
        c_list_y=[]
        errors=[]
        diffs=[]
        AC_calErrors=[]
        AL_calErrors=[]
        for o in collection:
            centroids,popt,perr, diff, calErr_AC,calErr_AL=o.evaluateCentroid_curveFit() # evaluates centroid of image
            errors.append(perr)
            c_list_x.append(centroids[0])
            c_list_y.append(centroids[1])
            diffs.append(diff)
            AC_calErrors.append(calErr_AC)
            AL_calErrors.append(calErr_AL)
        
        errors=np.array(errors)
        diffs=np.array(diffs)
        AC_calErrors=np.array(AC_calErrors)
        AL_calErrors=np.array(AL_calErrors)
        
        err_amp=np.mean(errors[:,0])
        err_x0=np.mean(errors[:,1])
        err_y0=np.mean(errors[:,2])
        err_sigmax=np.mean(errors[:,3])
        err_sigmay=np.mean(errors[:,4])
        err_theta=np.mean(errors[:,5])
        err_offset=np.mean(errors[:,6])
        
        max_x=np.max(c_list_x)
        min_x=np.min(c_list_x)
        mean_x=np.mean(c_list_x)
        std_x=np.std(c_list_x)
        
        max_y=np.max(c_list_y)
        min_y=np.min(c_list_y)
        mean_y=np.mean(c_list_y)
        std_y=np.std(c_list_y)
    
        print "Centroids estimated successfully"
        print "> min x0: " + str(min_x)
        print "> max x0: " + str(max_x)
        print "> mean x0: " + str(mean_x)
        print "> stdv x0: " + str(std_x)
        print
        print "> min y0: " + str(min_y)
        print "> max y0: " + str(max_y)
        print "> mean y0: " + str(mean_y)
        print "> stdv y0: " + str(std_y)
        print
        print "> mean error (pixel per pixel): " + str(round(np.mean(diffs)*100.0,3)) + "%"
        print "> mean calError (AC direction): " + str(np.mean(AC_calErrors)) + " micron"
        print "> mean calError (AL direction): " + str(np.mean(AL_calErrors)) + " micron"
        print
        print "> errors:"
        print ">> amplitude: " + str(err_amp) 
        print ">> x0: " +  str(err_x0) 
        print ">> y0: " + str(err_y0) 
        print ">> sigma_x: " + str(err_sigmax)
        print ">> sigma_y: " + str(err_sigmay)
        print ">> theta: " + str(err_theta)
        print ">> offset: "+ str(err_offset)
        print
    
    def shiftWindows(self, collection):
        dx=[]
        dy=[]
        values=[]
        for o in collection:
            rw=o.shiftWindow(o.centroidCF) # it returns the window with the relative coordinates wrt centroid
            # creating the (x,y,z) for the stacked version of the images (2D histogram)
            dx.append(rw[:,:,2])
            dy.append(rw[:,:,1])
            values.append(rw[:,:,0])
            
        # mapping the lists into numpy arrays   
        dx=np.array(dx)
        dy=np.array(dy)
        values=np.array(values)
        
        return dx,dy,values
    
    def divideTimeIntervals(self,hours):
        n_bins = int(self.time_interval/hours)+1
        bins = range(n_bins)
        collections = [[] for i in bins]
        for o in self.observations:
            for b in bins:
                if o.timestamp-self.first_obs < (b+1)*hours*60*60*1000*1000*1000:
                    o.timeBin=b
                    collections[b].append(o)
                    break
                    
        print "Created %d bins, each lasting %d hours" %(len(collections),hours)
        iterator=0
        for c in collections:
            print "> Bin %d has %d elements" %(iterator,len(c))
            iterator+=1
        return collections
    
            
    def evaluate_PCA(self, collection, from_library=True,normalized=False):
        
        from sklearn import decomposition
        
        input_data=[]
    
        for o in collection:
            flat_window=o.window.flatten()
            input_data.append(flat_window)
            
        input_data=np.matrix(input_data)
        
        if normalized==True:
            #normalizing data
            in_mean=np.mean(input_data,0)
            #in_stdv=np.std(input_data,0)
            data = (input_data-in_mean)
            #x_norm = (input_data-in_mean)/in_stdv
        
        else:
            data=input_data
            
        if from_library==True:
            # PCA using sklearn
            estimator=decomposition.PCA()
            estimator.fit(data)
            comp=estimator.components_
            
        else:
            # hand-made PCA
            x_t=np.transpose(input_data)
            cov=x_t*input_data
            ew,comp=np.linalg.eig(cov)
            
        return comp
            
    
        
            