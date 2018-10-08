'''
Created on 06 lug 2018

@author: jimmijamma
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import figure
from scipy import interpolate

S=4.223

if __name__ == '__main__':
    fp = open('scanratepar.dat','r')
    lines = fp.readlines()
    fp.close()
    t=[] # (days from J2000)
    sun_long=[] # (deg)
    sun_velocity=[] # (rad/sec)
    Omega=[] # (rad)
    d_omega=[] # (rad/sec)
    t_min=float(lines[1].split()[0])
    
    gamma=106.5 # Basic Angle (deg)
    
    t_list=[]
    zeta_list=[]
    eta_list=[]
    AC_smearing_list=[]
    AL_smearing_list=[]
    omega_list=[]
    for l in lines[1:]:
        values=np.array(l.split()).astype(float)
        # reading values
        t=(values[0]-t_min)*24
        t_list.append(t)
        sun_long=values[1]
        sun_velocity=values[2]
        Omega=values[3]
        omega_list.append(Omega)
        d_omega=values[4]
        
        _z=S*sun_velocity
        wz=60 # arcsec/sec
        zeta=np.radians(0.36) # (AC offset) = distance (deg) from the centre of the FOV
        AL_offset=0 # = distance (deg) from the centre of the FOV
        phi=np.radians(gamma*0.5+AL_offset) # celestial coordinate of the star
        
        _zeta=np.rad2deg(-_z*np.sin(Omega+phi))*3600000 # mas/sec
        _eta=np.rad2deg(_z*np.cos(Omega+phi)*np.tan(zeta))*3600000-wz # mas/sec
        
        zeta_list.append(_zeta)
        eta_list.append(_eta)
        
        integration_period=4.42 # (sec)
        pix_AC_size=176.8 # mas
        pix_AL_size=58.9 # mas
        AC_smearing=integration_period*abs(_zeta)/pix_AC_size # approximation of AC smearing (pix)
        AL_smearing=integration_period*abs(_eta+60)/pix_AL_size # approximation of AL smearing (pix)
        AC_smearing_list.append(AC_smearing)
        AL_smearing_list.append(AL_smearing) 
    
    
    new_xdata = np.linspace(min(t_list), max(t_list), len(t_list*1000))
    f = interpolate.interp1d(t_list, zeta_list, kind='cubic')
    new_ydata = f(new_xdata) 
    

    #plt.plot(new_xdata,new_ydata)
    #plt.savefig('kk.png')
    
    a=None
    b=None
    l_max=-1
    err_list=[]
    for ii,z in enumerate(new_ydata[:-4419]):
        lol=abs(abs(new_ydata[ii+4419])-(abs(z))) # difference of AC rate during integration period (mas/s)
        #err=lol*4.42/pix_AC_size
        #err_list.append(lol)
        if ii%3600000==0:
            print new_xdata[ii]
        if lol>l_max
            l_max=lol
            a=ii
            b=ii+4419
            
    print l_max, new_ydata[a],new_ydata[b]
    
    low=abs(new_ydata[a])*integration_period
    up=abs(new_ydata[b])*integration_period
    
    yes_integration=0
    for z in new_ydata[a:b+1]:
        yes_integration+=z
    yes_integration=yes_integration*0.001
    print yes_integration
    
    if a>3000:
        scarto=3000
    else: scarto=a
    
    ones=np.ones(len(new_ydata[a:b]))
    
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(131)
    ax1.plot(map(abs,new_xdata[a-scarto:a]),map(abs,new_ydata[a-scarto:a]),color='steelblue', lw=1)
    ax1.plot(new_xdata[a:b],new_ydata[a:b],color='green', lw=1.3)
    ax1.fill_between(map(abs,new_xdata[a:b]), 0, map(abs,new_ydata[a:b]), facecolor='green', alpha=0.5, label='smearing = %.3f mas'%yes_integration)
    ax1.plot(map(abs,new_xdata[b:b+scarto]),map(abs,new_ydata[b:b+scarto]),color='steelblue', lw=1)
    ax1.legend()
    ax1.set_title('With integration')
    ax2 = fig.add_subplot(132)
    ax2.plot(map(abs,new_xdata[a-scarto:a]),map(abs,new_ydata[a-scarto:a]),color='steelblue', lw=1)
    ax2.plot(new_xdata[a:b],new_ydata[a:b],color='green', lw=0.5)
    ax2.fill_between(map(abs,new_xdata[a:b]), 0, map(abs,new_ydata[a:b]), facecolor='green', alpha=0.2)
    ax2.plot(new_xdata[a:b],np.ones(len(new_ydata[a:b]))*new_ydata[a],color='indianred', lw=1.3)
    ax2.fill_between(map(abs,new_xdata[a:b]), 0, map(abs,np.ones(len(new_ydata[a:b]))*new_ydata[a]), facecolor='indianred', alpha=0.7, label='smearing = %.5f mas'%low)
    ax2.plot(map(abs,new_xdata[b:b+scarto]),map(abs,new_ydata[b:b+scarto]),color='steelblue', lw=1)
    ax2.legend()
    ax2.set_title('No integration (rounded down)')
    ax3 = fig.add_subplot(133)
    ax3.plot(map(abs,new_xdata[a-scarto:a]),map(abs,new_ydata[a-scarto:a]),color='steelblue', lw=1)
    ax3.plot(new_xdata[a:b],new_ydata[a:b],color='green', lw=0.5)
    ax3.fill_between(map(abs,new_xdata[a:b]), 0, map(abs,new_ydata[a:b]), facecolor='green', alpha=0.2)
    ax3.plot(new_xdata[a:b],np.ones(len(new_ydata[a:b]))*new_ydata[b],color='indianred', lw=1.3)
    ax3.fill_between(map(abs,new_xdata[a:b]), 0, map(abs,np.ones(len(new_ydata[a:b]))*new_ydata[b]), facecolor='indianred', alpha=0.7, label='smearing = %.3f mas'%up)
    ax3.plot(map(abs,new_xdata[b:b+scarto]),map(abs,new_ydata[b:b+scarto]),color='steelblue', lw=1)
    ax3.legend()
    ax3.set_title('No integration (rounded up)')
    fig.text(0.5, 0.04, 'time [h]', ha='center')
    fig.text(0.04, 0.5, '|AC rate| [mas/s]', va='center', rotation='vertical')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.savefig('situa.png')
    plt.close()
    
    dx=0.001
    dy = np.diff(map(abs,new_ydata))/dx
    
    plt.plot(new_xdata[:-1],dy)
    plt.savefig('derivative.png')
    
    '''
    # plot AC smearing
    w, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,w))
    axes = [ax, ax.twinx()]
    fig.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(new_xdata[:-4419],err_list, color='steelblue', lw=0.7)
    axes[0].tick_params(axis='y', colors='steelblue')
    axes[0].set_ylabel('AC smearing error [pix]')
    axes[1].plot(new_xdata[:-4419],new_ydata[:-4419], color='indianred', lw=0.7, ls=':')
    axes[1].tick_params(axis='y', colors='indianred')
    axes[1].set_ylabel('AC rate')
    plt.savefig('ACsmearingERROR.png')
    plt.close()
    '''
    
    '''
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
    
    # plot AC smearing
    w, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,w))
    axes = [ax, ax.twinx()]
    fig.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(t_list,AC_smearing_list, color='steelblue', lw=0.7)
    axes[0].tick_params(axis='y', colors='steelblue')
    axes[0].set_ylabel('AC smearing [pix]')
    axes[1].plot(t_list,zeta_list, color='indianred', lw=0.7, ls=':')
    axes[1].tick_params(axis='y', colors='indianred')
    axes[1].set_ylabel('AC rate')
    plt.savefig('ACsmearing.png')
    plt.close()
    
    # plot AL smearing
    w, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,w))
    axes = [ax, ax.twinx()]
    fig.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(t_list,AL_smearing_list, color='steelblue', lw=0.7)
    axes[0].tick_params(axis='y', colors='steelblue')
    axes[0].set_ylabel('AL smearing [pix]')
    axes[1].plot(t_list,eta_list, color='indianred', lw=0.7, ls=':')
    axes[1].tick_params(axis='y', colors='indianred')
    axes[1].set_ylabel('AL rate')
    plt.savefig('ALsmearing.png')
    plt.close()
    
    '''
        
        

    
    
    
        
    
        