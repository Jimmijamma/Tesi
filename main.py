'''
Created on 13 gen 2018

@author: jimmijamma
'''

import matplotlib as mpl 
mpl.use('TkAgg')
from Wrapper import Wrapper
from FrequencyAnalysis import FrequencyAnalysis


    
if __name__ == '__main__':
    
    print
    w=Wrapper('CalWrapper.json')
    fa=FrequencyAnalysis(w)
    collection=w.getCollection()
    
    fa.experiment_with_resize_PCA(collection,x_dim=180,y_dim=120)
    
    #fa.experiment_with_resize_aspect(collection, x_dim=1800, y_dim=1200)
    
    #fa.experiment(collection, x_dim=1800, y_dim=1200)
    
    fa.readResultsCooccurrence()
    fa.readResultsMoments()
    
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
    

    