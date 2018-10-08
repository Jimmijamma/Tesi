'''
Created on 15 lug 2018

@author: jimmijamma
'''
from subprocess import call,Popen,PIPE, STDOUT
import time
from Wrapper import Wrapper
from FrequencyAnalysis import FrequencyAnalysis
from ImgProc import ImgProc
import numpy as np
from matplotlib import pyplot as plt,figure
from scipy import interpolate



def buildFortranRoutine():
    call(['gfortran' ,'calcscanrate.for', '-o', 'calcscanrate.out'])
    time.sleep(2)
    

def callFortranRoutine(obmt_rev,dt,tmax):
    #print "Calling fortran routine for a period of %f days and dt of %f days" %(tmax,dt)
    p=Popen(['./calcscanrate.out'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    input_string=str(obmt_rev)+','+str(dt)+','+str(tmax)
    sys_out=p.communicate(input=input_string)[0]
    fp = open('scanratepar.dat','r')
    lines = fp.readlines()
    fp.close()
    return lines

def obmt2rev(obmt):
    one_rev=216e11
    obmt_rev=obmt*1.0/one_rev
    return obmt_rev

def milliseconds2days(ms):
    days=ms*1.0/86400000
    return days

def plotEtaZeta(t_list, eta_list, zeta_list):
    # plot AC rate
    w, h = figure.figaspect(6.)
    fig, axes = plt.subplots(figsize=(h,w))
    axes.set_xlabel('hours')
    axes.plot(t_list,zeta_list, color='steelblue', lw=0.7)
    axes.set_ylabel('AC rate [mas/s]')
    plt.savefig('zetaVStime.png')
    plt.close()
    
    # plot AL rate
    w, h = figure.figaspect(6.)
    fig, axes = plt.subplots(figsize=(h,w))
    axes.set_xlabel('hours')
    axes.plot(t_list,eta_list, color='steelblue', lw=0.7)
    axes.set_ylabel('AL rate [mas/s]')
    plt.savefig('etaVStime.png')
    plt.close()

def getImageSpeed(lines,FOV,ccsStrip,ccdRow):
    gamma=106.5 # Basic Angle (deg)
    S=4.223
    
    zeta_list=[]
    eta_list=[]
    t_list=[]
    omega_list=[]
    #t_min=float(lines[1].split()[0])

    for l in lines[1:]:
        values=np.array(l.split()).astype(float)
        #t=(values[0]-t_min)*24
        t=values[0]
        t_list.append(t)
        # reading values
        #sun_long=values[1] # (deg) 
        sun_velocity=values[2] # (rad/sec)
        Omega=values[3] # (rad)
        omega_list.append(Omega)
        d_omega=values[4] # (rad/sec)
        
        _z=S*sun_velocity
        wz=60000 # mas/sec
        zeta=np.radians(0.36) # (AC offset) = AC distance (deg) from the centre of the FOV
        eta=np.radians(0) # = AL distance (deg) from the centre of the FOV
        phi=np.radians(gamma*0.5+eta) # celestial coordinate of the star
        
        _zeta=np.rad2deg(-abs(_z)*np.sin(Omega+phi))*3600000 # mas/sec
        _eta=np.rad2deg(abs(_z)*np.cos(Omega+phi)*np.tan(zeta))*3600000-wz # mas/sec
        
        zeta_list.append(_zeta)
        eta_list.append(_eta)
        
    #plotEtaZeta(t_list, eta_list, zeta_list)
    
    return t_list, zeta_list, eta_list

def getSmearing(zeta_list, eta_list):
    dt=0.001 # 1 millisecond
    AC_integration=sum(map(abs,zeta_list))*dt
    AL_integration=sum(map(abs,[e+60 for e in eta_list]))*dt
    
    return AC_integration, AL_integration

def findInterval(t_list, zeta_list, start, stop):
    for ii,t in enumerate(t_list):
        if t>start:
            ii_start=ii
            break
    for ii,t in enumerate(t_list[ii_start:]):
        if t>stop:
            ii_stop=ii
            break
        
    
if __name__ == '__main__':
    #buildFortranRoutine()
    w=Wrapper('CalWrapper.json')
    ip=ImgProc()
    
    o=w.getCollection()[0]
    obmt=o.timestamp
    obmt_rev=obmt2rev(obmt)
    
    dt=milliseconds2days(1000)
    lines=callFortranRoutine(obmt_rev,dt=dt,tmax=1)
    print lines[1]
    #lines=out.splitlines()
    t_list, zeta_list, eta_list, omega_list=getImageSpeed(lines, w.fov, o.ccdStrip, o.ccdRow)
    
    wdt, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,wdt))
    ax.plot(t_list,map(abs,zeta_list),color='steelblue', lw=1)
    ax.set_ylabel('abs(zeta)')
    ax.set_xlabel('OMBT revolutions')
    ii=0
    for o in w.getCollection():
        if ii%100 == 0:
            hour=(o.timestamp-w.getCollection(1)[0].timestamp)*1.0/(1000*1000*1000*60*60)
            ot=o.integration_time
            obmt=o.timestamp
            obmt_rev=obmt2rev(obmt)
            
            dt=milliseconds2days(1000)
            tmax=milliseconds2days(ot)
            
            lines=callFortranRoutine(obmt_rev,dt=dt,tmax=0.01)
            print lines[1]

            #lines=out.splitlines()
            _t, _z, _e, _o=getImageSpeed(lines, w.fov, o.ccdStrip, o.ccdRow)
            start=min(_t)
            stop=max(_t)
            ax.plot(_t,map(abs,_z),color='green', lw=1.3)
            ax.fill_between(_t,map(abs,_z), facecolor='green', alpha=0.5)
        ii+=1
    plt.savefig('ACrateIntervals.png')
    plt.close()
    