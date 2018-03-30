'''
Created on 13 gen 2018

@author: jimmijamma
'''
import json
import numpy as np
import matplotlib as mpl 
from scipy.linalg import decomp
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Wrapper import Wrapper,pix_size_AC,pix_size_AL
import png
from scipy import optimize as opt
from skimage.feature import register_translation
from sklearn import decomposition

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
    
    input_matrix=[]
    
    for o in w.observations:
        flat_window=o.window.flatten()
        input_matrix.append(flat_window)
        
    input_matrix=np.matrix(input_matrix)
    
    (N, F)=np.shape(input_matrix)
    print N,F
    
    
    
    ''''
    #normalizing data

    in_mean=np.mean(input_matrix,0)
    in_stdv=np.std(input_matrix,0)

    x_norm = (input_matrix-in_mean)
    
    # hand-made
    x_t=np.transpose(input_matrix)
    cov=x_t*input_matrix
    ew,comp=np.linalg.eig(cov)
    
    
    cumsum=sum(ew)
    l_sum=0
    for l in range(len(ew)):
        l_sum+=ew[l]
        if l_sum>0.999*cumsum:
            break
    print "dimension: " + str(l)
    '''
    
    # PCA using sklearn
    estimator=decomposition.PCA()
    estimator.fit(input_matrix)
    comp=estimator.components_
    
    for i in range(8):
        eigenimage=(np.reshape(comp[i,:], (12,18), order=0)-1)*65000/2+65000
        eigenimage=np.array(eigenimage)
        shift, error, diffphase = register_translation(eigenimage,w.observations[0].window,1000)
        print shift[0]+5.5
        print shift[1]+8.5
        print [w.observations[0].calCentroid_AC,w.observations[0].calCentroid_AL]
        #if i==1:
        #    eig_sample=eigenimage
        #    print eig_sample
        f=open('eigenimages/eigenimage'+str(i)+'.png', 'wb')
        win = map(np.uint16,eigenimage)
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

    
    '''
    h=np.nan_to_num(h)
    xe=xe[:-1]
    ye=ye[:-1]
    xe, ye = np.meshgrid(xe, ye)
        
    # initial guess of parameters (amplitude,dx,dy,sigma_x,sigma_y,theta,offset)
    initial_guess = (5000,0,0,1,1,0,1)
    # constraints for fitting ((lower bounds),(upper bounds))
    bounds=((0,-9,-6,-np.inf,-np.inf,0,-np.inf),(500000,9,6,np.inf,np.inf,np.inf,np.inf))
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (xe, ye), h.ravel(),p0=initial_guess, bounds=bounds,maxfev=1000000)
    
    g=twoD_Gaussian((xe, ye), popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(ye,xe,g.reshape(119,179)) 
    fig.savefig('bhoooo.png')
    
    '''
    
    
 
    
    
    