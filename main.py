'''
Created on 13 gen 2018

@author: jimmijamma
'''

import matplotlib as mpl 
mpl.use('TkAgg')
from Wrapper import Wrapper
from FrequencyAnalysis import FrequencyAnalysis
from ImgProc import ImgProc
import numpy as np
import cv2
from matplotlib import pyplot as plt
from matplotlib import figure
from scipy.optimize import curve_fit
from progress.bar import ChargingBar, Bar
import png

Z=[ 1.04477081, -1.57731015, -0.13579937, 1.10002968]

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
if __name__ == '__main__':
    
    
    print
    w=Wrapper('CalWrapper.json')
    fa=FrequencyAnalysis(w)
    ip=fa.improc
    '''
    img=w.observations[0].window
    img=ip.img_interpolateImage(img, 180, 120)
    img = (img/256)
    f=open('original.png', 'wb')
    win=np.array(img)
    cv2.imwrite('original.png',win)
    
    img = cv2.imread('original.png')
    res=cv2.imencode('.png',img)
    cv2.imwrite('original.png',img)
    img=res[1]
    plt.savefig('hist.png')
    ret, thresh = cv2.threshold(img,200 , 255, 0)
    l_img=[]
    for el in thresh:
        l_img.append(el[0])
    print l_img
    cv2.imwrite('thresh.png',thresh)

    im2, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    print contours
   '''
    
    
    
    # z=fa.polyfit_ACrate_ROIaspect(w.observations, deg=3, x_dim=180, y_dim=120)
    z=Z
    
    fa=FrequencyAnalysis(w)
    ip=fa.improc
    collection=w.observations
    fa.experiment_with_polyfit(collection, coeff=Z, x_dim=720, y_dim=480)
    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    
    
    fa=FrequencyAnalysis(w)
    ip=fa.improc
    collection=w.observations
    fa.experiment_with_resize_aspect(collection, x_dim=720, y_dim=480)
    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    
    fa=FrequencyAnalysis(w)
    ip=fa.improc
    collection=w.observations
    fa.experiment_with_resize_PCA(collection, x_dim=180, y_dim=120)
    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    '''      
    fa.display_timedomain(l_hours, lol, 'ROI_aspect', l_ACrate, dir_name='.')
    '''        
    '''
    fig,ax =plt.subplots()
    axes=[ax,ax.twinx()]
    axes[0].plot(lol, 'go', ms=0.6)
    axes[0].tick_params(axis='y', colors='green')
    axes[1].plot(l_ACrate, color='red')
    plt.savefig('miao.png')
    plt.close()
    '''
       
    ''' 
    fa=FrequencyAnalysis(w)
        
    collection=w.getCollection()
    
    fa.experiment_with_resize_PCA(collection,x_dim=360,y_dim=240)
    
    #fa.experiment_with_resize_aspect(collection, x_dim=1800, y_dim=1200)
    
    #fa.experiment(collection, x_dim=1800, y_dim=1200)
    
    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    '''
    
    '''
    w.divideTimeIntervals(3)
    
    # create images of the set of observations
    for o in w.observations:
        o.createImage(folder, filename, img_format) 
    

    # evaluating PCA
    n_components=10
    from_library=False
    comp=w.evaluate_PCA(w.getCollection(), from_library=True, normalized=False, n_components=n_components)
    
    
    
    cllctn=[]
    for o in w.observations:
        cllctn.append(o.window)
    
    mean_img= np.mean(cllctn, axis=0)

    
    f=open('mean_img.png', 'wb')
    win = map(np.uint16,mean_img)
    writer = png.Writer(width=len(win[0]), height=len(win), bitdepth=16, greyscale=True)
    writer.write(f, win)
    f.close()
    
    new_collection=w.filterBadRegistration(w.observations)
    
    
    
    rel_shift=[]
    coor=[0,0]
    for o in new_collection:
        if o.id>0:
            rel_shift.append([o.calCentroid_AC-coor[0],o.calCentroid_AL-coor[1]])
        else:
            rel_shift.append(coor)
        coor=[o.calCentroid_AC,o.calCentroid_AL]
        
    
    rel_shift_B=[]
    lol=[0,0]
    for i in range(len(new_collection)):
        
        shift, error, diffphase = register_translation(mean_img,new_collection[i].window,1000)
        if i>0:
            rel_shift_B.append([shift[0]-lol[0],shift[1]-lol[1]])
        else:
            rel_shift_B.append(lol)
        lol=shift
            
    
    errs_AC=[]
    errs_AL=[]
    for i in range(1,len(new_collection)):
        errs_AC.append(abs(rel_shift_B[i][0]-rel_shift[i][0]))
        errs_AL.append(abs(rel_shift_B[i][1]-rel_shift[i][1]))
        
    print my_util.AC_pix2micron((np.mean(errs_AC)))
    print my_util.AL_pix2micron((np.mean(errs_AL)))
    print
    
    
    # estimating centroids for each transit fitting with a Gaussian curve
    w.estimateCentroids_fit()
    
    # shifting the windows wrt their centroid
    dx,dy,values=w.shiftWindows()
    
    # creating a 2D histogram
    x_res=180 # resolution on x-axis
    y_res=120
    xe,ye,h=w.create2Dhistogram('2Dhist_120x180.png', dx, dy, values, x_res, y_res, normed=True)
 
    '''
    

    