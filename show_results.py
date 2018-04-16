'''
Created on 09 apr 2018

@author: jimmijamma
'''
from matplotlib import pyplot as plt
from matplotlib import figure
import numpy as np
from scipy import interpolate

class Measure(object):
    def __init__(self,id,hours,area_roi,aspect_roi,R,third_moment,uniformity):
        self.id=id
        self.hours=hours
        self.area_roi=area_roi
        self.aspect_roi=aspect_roi
        self.R=R
        self.third_moment=third_moment
        self.uniformity=uniformity
        
def interpolate_measurements(xdata,ydata):
    new_xdata = np.linspace(min(xdata), max(xdata), len(xdata))
    f = interpolate.interp1d(xdata, ydata)
    new_ydata = f(new_xdata)   # use interpolation function returned by `interp1d`
    return new_xdata,new_ydata
    
def evaluate_fft(x_data,y_data,y_name,peak_threshold):
    
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
    
    display_fft(periods,pos_spectrum, param_name=y_name)


    return pos_freq,pos_spectrum


def display_fft(periods,spectrum,param_name,peaks=None):
    
    fig,ax=plt.subplots()
    ax.plot(periods,spectrum, color='green')
    ax.set_xlabel('Period [h]')
    ax.set_ylabel('Amplitude')
    plt.ylim([None, max(spectrum)*1.1])
    '''
    bbox_props = dict(boxstyle='round', alpha=0.7, ec="b", lw=1, fc='white')
    max_x=np.max(periods)
    max_y=np.max(spectrum)
    for i,s in enumerate(spectrum):
        if s>=0.05*max(spectrum):
            ax.text(periods[i], spectrum[i], "T="+str(round(periods[i],3)), ha="center", va="center",
            size=10, bbox=bbox_props)
    '''
    plt.title("Periodicity of '%s' parameter" % param_name)
    plt.savefig('fft_'+param_name+'.png')
    
def display_timedomain(x_var,y_var,y_name):
    
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
    
    plt.savefig('ACmotionVS'+y_name+'.png')
    
        
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
    
    
    vars_dict=dict({'area_roi':area_roi,'aspect_roi':aspect_roi,'R':R,'third_moment':third_moment,'uniformity':uniformity,})
    for v in vars_dict:
        xvar=hours
        yvar=vars_dict[v]
        interp_x,interp_y=interpolate_measurements(xvar, yvar-np.mean(yvar))
        display_timedomain(interp_x, interp_y, v)
        freq,spectrum=evaluate_fft(interp_x, interp_y, v,peak_threshold=0.05)
        
        
 
    
        