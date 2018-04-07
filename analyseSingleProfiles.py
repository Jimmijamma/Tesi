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
import matplotlib.pyplot as plt

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
    
    ref_images=[]
    # create images of the set of observations
    for o in w.observations:
        o.createImage(folder, filename, img_format) 
        ref_images.append(o.imgpath)
    
    # Limit starting size
    images_orig = [cv2.resize(imageio.imread(im), (18,12)) for im in ref_images[:3]] 
    # interpolation methods to compare
    methods=[("area", cv2.INTER_AREA), 
         ("nearest", cv2.INTER_NEAREST), 
         ("linear", cv2.INTER_LINEAR), 
         ("cubic", cv2.INTER_CUBIC), 
         ("lanczos4", cv2.INTER_LANCZOS4)]
    
    image_set = [[cv2.resize(im, (1800,1200), interpolation=m[1]) for m in methods] for im in images_orig]
    print np.shape(image_set)
    display(image_set)
        
    
    
    
        
    
    
    