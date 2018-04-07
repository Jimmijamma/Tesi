'''
Created on 07 apr 2018

@author: jimmijamma
'''

from MyUtil import MyUtil
import Wrapper
import json
from os import listdir,makedirs
from os.path import isfile,join


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
    
    # create images of the set of observations
    for o in w.observations:
        o.createImage(folder, filename, img_format) 
        
    for o in w.observations:
        print o.imgpath
    
        
    
    
    