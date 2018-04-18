'''
Created on 07 apr 2018

@author: jimmijamma
'''

from MyUtil import MyUtil
import Wrapper
import json
import cv2
import imageio
import numpy as np
import itertools
import matplotlib.gridspec as gridspec
import pandas as pd
import time
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log, blob_doh, greycomatrix, greycoprops
from skimage.color import rgb2gray

import matplotlib.pyplot as plt
from dask.dataframe.tests.test_rolling import idx
from conda._vendor.toolz.itertoolz import second

# function to display images
def display(images, titles=['']):
    if isinstance(images[0], list):
        c = len(images[0])
        r = len(images)
        images = list(itertools.chain(*images))
    else:
        c = len(images)
        r = 1
    fig=plt.figure(figsize=(4*c, 4*r))
    gs1 = gridspec.GridSpec(r, c, wspace=0, hspace=0)
    #gs1.update(wspace=0.01, hspace=0.01) # set the spacing between axes. 
    titles = itertools.cycle(titles)
    for i in range(r*c):
        im = images[i]
        title = titles.next()
        plt.subplot(gs1[i])
        # Don't let imshow doe any interpolation
        plt.imshow(im, cmap='gray', interpolation='none')
        plt.axis('off')
        if i < c:
            plt.title(title)
    plt.tight_layout()
    fig.savefig('InterpolationComparison.png')
    
def detectROI(image,zero_padding=True):
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

def analyseROI(ROI):
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
    
    
    '''
    print "R: %f" % R
    print "third_moment: %f" % third_moment
    print 'uniformity: %f' % uniformity
    print 'entropy: %f' % entropy
    '''
    
    
    
    
    

    

if __name__ == '__main__':
    my_util=MyUtil()
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
    w=Wrapper.Wrapper(json_obj)
    
    '''
    ref_images=[]
    # create images of the set of observations
    for o in w.observations:
        o.createImage(folder, filename, img_format) 
        ref_images.append(o.imgpath)
        

    
    # Limit starting size
    images_orig = [imageio.imread(im) for im in ref_images[:3]] 
    # interpolation methods to compare
    methods=[("area", cv2.INTER_AREA), 
         ("nearest", cv2.INTER_NEAREST), 
         ("linear", cv2.INTER_LINEAR), 
         ("cubic", cv2.INTER_CUBIC), 
         ("lanczos4", cv2.INTER_LANCZOS4)]
    
    image_set = [[cv2.resize(im, (1800,1200), interpolation=m[1]) for m in methods] for im in images_orig]
    display(image_set)
    
    
    # Read image
    img_orig = cv2.imread(ref_images[1], cv2.IMREAD_ANYCOLOR)
    im=cv2.resize(img_orig,(1800,1200),interpolation=cv2.INTER_LANCZOS4)
    # Set up the detector with default parameters.
    '''
    
    collection=w.observations
    txt_file = open('results.txt','w') 
    for o in collection:     
        # original image
        image=o.window
        # interpolated image
        interp_image=cv2.resize(image,(1800,1200),interpolation=cv2.INTER_LANCZOS4)
        # detecting blob
        ROI, row_interval, col_interval=detectROI(interp_image)
        
        glcm = greycomatrix(ROI.astype(int)/256, [1], [0], 256, symmetric=True, normed=True)
        contrast = greycoprops(glcm, prop='contrast')[0, 0]
        dissimilarity = greycoprops(glcm, prop='dissimilarity')[0, 0]
        homogeneity = greycoprops(glcm, prop='homogeneity')[0, 0]
        energy = greycoprops(glcm, prop='energy')[0, 0]
        correlation = greycoprops(glcm, prop='correlation')[0, 0]
        ASM = greycoprops(glcm, prop='ASM')[0, 0]
        
        #print [contrast, dissimilarity, homogeneity, energy, correlation, ASM]
        
        hours=(o.timestamp-collection[0].timestamp)*1.0/(60*60*1000*1000*1000)
        
        #area_roi,aspect_roi,R,third_moment,uniformity=analyseROI(ROI)
        
        
        
        #txt_file.write('%d %.6f %d %d %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACmotion,area_roi,aspect_roi,R,third_moment,uniformity))
        txt_file.write('%d %.6f %d %6f %.6f %.6f %.6f %.6f %.6f\n' % (o.id,hours,o.ACmotion,contrast, dissimilarity, homogeneity, energy, correlation, ASM))
        
        print "%d of %d" % (o.id,len(collection))
        '''
        # displaying things
        fig, ax = plt.subplots(3,1)
        plt.subplots_adjust(hspace=0.6)
        ax[0].set_title('Original image (12x18)')
        ax[0].imshow(image,cmap='gray')
        ax[1].set_title('Interpolated image (1200x1800)')
        ax[1].imshow(interp_image,cmap='gray')
        ax[2].set_title('Blob detection')
        ax[2].imshow(interp_image,cmap='gray')
        ax[2].add_patch(patches.Rectangle((col_interval[0],row_interval[0]),np.shape(ROI)[1],np.shape(ROI)[0], fill=False, edgecolor='red'))
        fig.savefig('triplot.png')
        '''
    txt_file.close()
    
    
    '''
    prms = cv2.SimpleBlobDetector_Params()
    prms.filterByColor=True
    prms.blobColor=255
    detector = cv2.SimpleBlobDetector_create(prms)
     
    # Detect blobs.
    keypoints = detector.detect(im)
    print keypoints[0].pt # coordinate of keypoint
    print keypoints[0].size # radius of keypoint
     
    # Draw detected blobs as red circles.
    # cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS ensures the size of the circle corresponds to the size of blob
    im_with_keypoints = cv2.drawKeypoints(im, keypoints, np.array([]), (0,0,255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
     
    # Show keypoints
    params=list()
    params.append(cv2.IMWRITE_PNG_COMPRESSION)
    params.append(0)
    cv2.imwrite('blob_img.png', im_with_keypoints, params)
    '''
    
    '''
    image_gray=im
    blobs_log = blob_log(image_gray, max_sigma=100, num_sigma=10, threshold=.05)
    
    # Compute radii in the 3rd column.
    blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
    
    blobs_dog = blob_dog(image_gray, max_sigma=100, threshold=.05)
    blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)
    
    blobs_doh = blob_doh(image_gray, max_sigma=100, threshold=.005)
    
    blobs_list = [blobs_log, blobs_dog, blobs_doh]
    colors = ['yellow', 'lime', 'red']
    titles = ['Laplacian of Gaussian', 'Difference of Gaussian',
              'Determinant of Hessian']
    sequence = zip(blobs_list, colors, titles)
    
    fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
    ax = axes.ravel()
    
    for idx, (blobs, color, title) in enumerate(sequence):
        ax[idx].set_title(title)
        ax[idx].imshow(im)
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
            ax[idx].add_patch(c)
        ax[idx].set_axis_off()

    fig.savefig('scikit.png')
    '''