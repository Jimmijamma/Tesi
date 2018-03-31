'''
Created on 13 gen 2018

@author: jimmijamma
'''
import json
import numpy as np
import matplotlib as mpl 
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Wrapper import Wrapper
import png
from scipy import optimize as opt
from skimage.feature import register_translation
from Observation import Observation,twoD_Gaussian


if __name__ == '__main__':
    # path and format for the images
    folder='/Users/jimmijamma/Desktop/bho'
    filename=None
    img_format='png'
    
    # loading the wrapper from JSON
    json_data=open('CalWrapper.json')
    json_obj = json.load(json_data)
    json_data.close()
    
    print
    # loading the data structure and calling methods
    w=Wrapper(json_obj)
    
    '''
    w.divideTimeIntervals(3)
    
    
    # create images of the set of observations
    for o in w.observations:
        o.createImage(folder, filename, img_format) 
    '''
    
    n_components=10
    
    from_library=False
    
    comp=w.evaluate_PCA(w.observations, from_library=False, normalized=False)
    
    '''
    for i in range(n_components):
        if from_library==True:
            eigenimage=(np.reshape(comp[i,:], (12,18), order=0)-1)*65000/2+65000
        else:
            eigenimage=(np.reshape(comp[:,i], (12,18), order=0)-1)*65000/2+65000
        eigenimage=np.array(eigenimage)
        shift, error, diffphase = register_translation(eigenimage,w.observations[0].window,1000)
        print shift[0]+5.5
        print shift[1]+8.5
        print [w.observations[i].calCentroid_AC,w.observations[i].calCentroid_AL]
        #if i==1:
        #    eig_sample=eigenimage
        #    print eig_sample
        f=open('eigenimages/eigenimage'+str(i)+'.png', 'wb')
        win = map(np.uint16,eigenimage)
        writer = png.Writer(width=len(win[0]), height=len(win), bitdepth=16, greyscale=True)
        writer.write(f, win)
        f.close()
    '''
    
    cllctn=[]
    for o in w.observations:
        cllctn.append(o.window)
    
    mean_img= np.mean(cllctn, axis=0)

    
    f=open('mean_img.png', 'wb')
    win = map(np.uint16,mean_img)
    writer = png.Writer(width=len(win[0]), height=len(win), bitdepth=16, greyscale=True)
    writer.write(f, win)
    f.close()
    
    '''
    err_AC=[]
    err_AL=[]
    for i in range(1,w.n_transits):
        result1=[w.observations[0].calCentroid_AC, w.observations[0].calCentroid_AL]
        result2=[w.observations[i].calCentroid_AC, w.observations[i].calCentroid_AL]
        calResult= [ -(result2[0]-result1[0]), -(result2[1]-result1[1])]
        shift, error, diffphase = register_translation(w.observations[0].window,w.observations[i].window,1000)
        if shift[0]<=1.0 and shift[1]<=1.0: 
            err_AC.append(abs(calResult[0]-shift[0]))
            err_AL.append(abs(calResult[1]-shift[1]))
        else:
            print "error " +str(i) + ": " + str(shift[0]) + " " + str(shift[1])
        print "AC: " + str(np.mean(err_AC)*pix_size_AC)
        print "AL: " + str(np.mean(err_AL)*pix_size_AL)
       
    '''
         
    '''
    # estimating centroids for each transit fitting with a Gaussian curve
    w.estimateCentroids_fit()
    
    # shifting the windows wrt their centroid
    dx,dy,values=w.shiftWindows()
    
    # creating a 2D histogram
    x_res=180 # resolution on x-axis
    y_res=120
    xe,ye,h=w.create2Dhistogram('2Dhist_120x180.png', dx, dy, values, x_res, y_res, normed=True)
 
    '''

    

    
    
 
    
    
    