'''
Created on 13 gen 2018

@author: jimmijamma
'''

import numpy as np
import png
from scipy import optimize as opt
from astropy.version import timestamp

pix_size_AL=47260/4500
pix_size_AC=60000/1966 # micron

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

class Observation(object):
    '''
    classdocs
    '''

    def __init__(self, wrapper, id, window, gating, transitid, timestamp, ACmotion, ACrate, calCentroid_AC=None, calCentroid_AL=None):
        '''
        Constructor
        '''
        self.wrapper=wrapper
        self.id=id
        self.calCentroid_AC=calCentroid_AC
        self.calCentroid_AL=calCentroid_AL
        self.window=window
        self.ccdRow=wrapper.ccdRow
        self.ccdStrip=wrapper.ccdStrip
        self.mag=wrapper.mag
        self.wavn=wrapper.wavn
        self.gating=gating
        self.timestamp=timestamp # ns
        self.transitid=transitid
        self.ACmotion=ACmotion
        if ACrate==None:
            ACrate=-1
        self.ACrate=ACrate
        self.imgpath=None
        self.centroidCM=None
        self.centroidCF=None
        self.totIntensity=None
        self.timeBin=None
        self.ROIaspect=None
                
    def createImage(self, folder, filename, format):
        if filename==None:
            filename= str(self.timestamp)
            total_path=folder + "/" + filename + "." + format
        else:
            total_path=folder + "/" + filename + str(self.id) + "." + format
            
        f=open(total_path, 'wb')
        win = map(np.uint16,self.window)
        writer = png.Writer(width=len(win[0]), height=len(win), bitdepth=16, greyscale=True)
        writer.write(f, win)
        f.close()
        self.imgpath=total_path
        
    # function that evaluates centroid by fitting the window with a curve
    def evaluateCentroid_curveFit(self):

        x = np.linspace(start=0, stop=17, num=18)
        y = np.linspace(start=0, stop=11, num=12)
        x, y = np.meshgrid(x, y)
        
        # initial guess of parameters (amplitude,dx,dy,sigma_x,sigma_y,theta,offset)
        initial_guess = (1500,8,5,1,1,0,1)
        # constraints for fitting ((lower bounds),(upper bounds))
        bounds=((0,0,0,-np.inf,-np.inf,0,-np.inf),(500000,17,11,np.inf,np.inf,np.inf,np.inf))
        
        popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), self.window.ravel(),p0=initial_guess, bounds=bounds,maxfev=1000000)
        
        # compute fitted curve
        g=twoD_Gaussian((x, y), popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
        
        # relative difference between original window and fitted one (absolute value)
        diffs=np.mean(np.abs(self.window.ravel()-g)/self.window.ravel())
        calErr_AC=abs(self.calCentroid_AC-popt[2])*pix_size_AC
        calErr_AL=abs(self.calCentroid_AL-popt[1])*pix_size_AL

        perr = np.sqrt(np.diag(pcov))
        self.centroidCF=[popt[1],popt[2]]
        return self.centroidCF,popt,perr,diffs,calErr_AC,calErr_AL
    
        
    # function that evaluates centroid by computing center of mass
    def evaluateCentroid_centerOfMass(self):
        b00=0
        # computing total intensity of image
        for r in self.window:
            b00=b00+sum(r)
        
        # computing the 1st order moments of image
        b10=0
        b01=0
        for y in range(len(self.window)):
            for x in range(len(self.window[0])):
                b10=b10+x*self.window[y][x]
                b01=b01+y*self.window[y][x]
                
        b10=b10/b00 # x centroid
        b01=b01/b00 # y centroid
        
        centroid=[b10,b01] # (x,y) centroid
        self.centroidCM=centroid
        self.totIntensity=b00
        return self.centroidCM
    
    
    # function that shifts the original window wrt its estimated centroid
    def shiftWindow(self, centroid):
        cx=centroid[0]
        cy=centroid[1]
        rel_window=[]
        for y in range(len(self.window)):
            row=[]
            for x in range(len(self.window[0])):
                dx=cx-x
                dy=cy-y
                row.append([self.window[y][x],dy,dx])
            rel_window.append(row)
            
        rel_window=np.array(rel_window)
        return rel_window