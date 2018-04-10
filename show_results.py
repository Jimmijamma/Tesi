'''
Created on 09 apr 2018

@author: jimmijamma
'''
from matplotlib import pyplot as plt
from matplotlib import figure
import numpy as np
from sympy.plotting.intervalmath.interval_arithmetic import interval

class Measure(object):
    def __init__(self,id,hours,area_roi,aspect_roi,R,third_moment,uniformity):
        self.id=id
        self.hours=hours
        self.area_roi=area_roi
        self.aspect_roi=aspect_roi
        self.R=R
        self.third_moment=third_moment
        self.uniformity=uniformity
        
        
if __name__ == '__main__':
    f=open('results.txt','r')
    measure_list=[]
    for line in f:
        l=line.split()
        m=[l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]]
        measure_list.append(m)
    f.close()
    
    ids=[int(row[0]) for row in measure_list]
    hours=[float(row[1]) for row in measure_list]
    ACmotion=[int(row[2]) for row in measure_list]
    area_roi=[int(row[3]) for row in measure_list]
    aspect_roi=[float(row[4]) for row in measure_list]
    R=[float(row[5]) for row in measure_list]
    third_moment=[float(row[6]) for row in measure_list]
    uniformity=[float(row[7]) for row in measure_list]
    
    '''
    plt.figure()
    plt.plot(hours,aspect_roi)
    plt.savefig('aspectvshours.png')
    '''
    w, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,w))

    # Twin the x-axis twice to make independent y-axes.
    axes = [ax, ax.twinx()]
    
    # Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.75)
     
    axes[0].set_xlabel('hours')
    axes[0].plot(hours,area_roi, color='green')
    axes[0].tick_params(axis='y', colors='green')
    axes[0].set_ylabel('ROI area [pix$^2$]', color='green')
    
    axes[1].plot(hours,ACmotion, color='red')
    axes[1].tick_params(axis='y', colors='red')
    axes[1].set_ylabel('AC motion',color='red')
    
    plt.savefig('Rvshours.png')
    
    
    A=np.fft.fft(area_roi)
    freq = np.fft.fftfreq(len(A))
    spectrum=np.abs(A)
    
    fig,ax=plt.subplots()
    ax.plot(freq[:len(freq)/2],spectrum[:len(spectrum)/2])
    plt.savefig('fft.png')
    threshold = 0.01 * max(spectrum)
    mask = abs(spectrum) > threshold
    
    
    start=0
    intervals=[]
    for stop in hours:
        intervals.append(stop-start)
        start=stop
    fs= 1.0/np.mean(intervals)
    
    peaks = freq[mask]*fs
    print peaks
        
    
 
    
        