'''
Created on 26 giu 2018

@author: jimmijamma
'''

import numpy as np

if __name__ == '__main__':
    q01,q02,q03,q04=[0.8,0.7,0.8,-0.5]
    q11,q12,q13,q14=[0.8,0.2,0.8,0.3]
    
    wx = 2*(q04*q11+q03*q12-q02*q13-q01*q14)
    wy = 2*(-q03*q11+q04*q12+q01*q13-q02*q14)
    wz = 2*(q02*q11-q01*q12+q04*q13-q03*q14)
    print (wx, wy, wz)
    
    gamma=np.radians(106.5)
    fov=-1
    
    zeta=fov*wx*np.sin(gamma/2)+wy*np.cos(gamma/2)
    print zeta