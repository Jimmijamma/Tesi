'''
Created on 15 ago 2018

@author: jimmijamma
'''
from Wrapper import Wrapper
from FrequencyAnalysis import FrequencyAnalysis
import cv2
import png
import numpy as np

if __name__ == '__main__':
    w=Wrapper('CalWrapper.json')
    fa=FrequencyAnalysis(w)
    ip=fa.improc
    
    '''
    im=w.observations[2222].window
    dict_cv={'linear':cv2.INTER_LINEAR, 'nearest':cv2.INTER_NEAREST, 'area':cv2.INTER_AREA, 'cubic':cv2.INTER_CUBIC, 'lanczos':cv2.INTER_LANCZOS4}
    for k in dict_cv:
        new_im=ip.img_interpolateImage(im, 1800, 1200, algorithm=dict_cv[k])
        print k
        
        f=open(str(k)+'.png', 'wb')
        win = map(np.uint16,new_im)
        writer = png.Writer(width=len(win[0]), height=len(win), bitdepth=16, greyscale=True)
        writer.write(f, win)
        f.close()
    '''
    
    #fa.experiment(w.observations, 180, 120)
    fa.experiment_with_resize_PCA(w.observations, 360, 240)
    #fa.experiment_with_resize_aspect(w.observations, x_dim=1800, y_dim=1200)
    fa.readResultsCooccurrence()
    fa.readResultsMoments()