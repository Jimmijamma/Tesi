'''
Created on 30 mar 2018

@author: jimmijamma
'''

class MyUtil(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        
        self.sensor_size_AL=47260 # micron
        self.n_pixel_AL=4500
        
        self.sensor_size_AC=60000 # micron
        self.n_pixel_AC=1966 
        
        
        
        def AL_micron2pix(micrometers):
            pix=micrometers*self.n_pixel_AL
            pix=pix/self.sensor_size_AL
            return pix
        
        def AC_micron2pix(micrometers):
            pix=micrometers*self.n_pixel_AC
            pix=pix/self.sensor_size_AC
            return pix
        
        def AL_pix2micron(n_pixels):
            micron=n_pixels*self.sensor_size_AL
            micron=micron/self.n_pixel_AL
            return micron
        
        def AC_pix2micron(n_pixels):
            micron=n_pixels*self.sensor_size_AC
            micron=micron/self.n_pixel_AC
            return micron
        
        
            