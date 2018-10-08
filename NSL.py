'''
Created on 22 ago 2018

@author: jimmijamma
'''
from subprocess import call,Popen,PIPE, STDOUT
import time
from Wrapper import Wrapper
from FrequencyAnalysis import FrequencyAnalysis
from ImgProc import ImgProc
import numpy as np
from matplotlib import pyplot as plt,figure
from scipy import interpolate,fft,ifft,conj,signal
import matplotlib

def periodic_corr(x, y):
    """Periodic correlation, implemented using the FFT.

    x and y must be real sequences with the same length.
    """
    return ifft(fft(x) * fft(y).conj()).real

def buildFortranRoutine_init():
    call(['gfortran' ,'calcscanwow.for', '-o', 'calcscanwow.out'])
    time.sleep(2)
    
def buildFortranRoutine():
    call(['gfortran' ,'calcscanrate.for', '-o', 'calcscanrate.out'])
    time.sleep(2)
    
def callFortranRoutine_init(dt,tmax):
    #print "Calling fortran routine for a period of %f days and dt of %f days" %(tmax,dt)
    p=Popen(['./calcscanwow.out'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    input_string='2000,1,1,0.0,'+str(dt)+','+str(tmax)
    print input_string
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

def getInitParams(w):
    o=w.getCollection()[0]
    #obmt=o.timestamp-(0.40*60*60*1000*1000*1000)
    obmt=o.timestamp
    #calling the fortran routine to get the initial date of the collection
    obmt_rev=obmt2rev(obmt)
    l=callFortranRoutine(obmt_rev, dt=0.1, tmax=1, Omega0=0, nu0=0)
    start_date=l[1].split()[0]
    print start_date
    
    #calling again the routine to get the initial Omega0 and nu0
    lines=callFortranRoutine_init(dt=start_date,tmax=start_date)
    init = lines[-1].split()
    print init
    Omega0=init[3]
    nu0=init[5]
    return obmt,Omega0,nu0

def callFortranRoutine(obmt_rev,dt,tmax,Omega0,nu0):
    #print "Calling fortran routine for a period of %f days and dt of %f days" %(tmax,dt)
    p=Popen(['./calcscanrate.out'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    input_string=str(obmt_rev)+','+str(dt)+','+str(tmax)+','+str(nu0)+','+str(Omega0)
    sys_out=p.communicate(input=input_string)[0]
    #lines=sys_out.splitlines()
    fp = open('scanratepar.dat','r')
    lines = fp.readlines()
    fp.close()
    return lines

def getImageSpeed(lines,FOV,ccsStrip,ccdRow):
    gamma=106.5 # Basic Angle (deg)
    S=4.223
    if FOV==1:
        sn=1
    else:
        sn=-1
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
        zeta=0 # (AC offset) = AC distance (deg) from the centre of the FOV
        eta=0 # = AL distance (deg) from the centre of the FOV
        phi=np.radians(sn*gamma*0.5+eta) # celestial coordinate of the star
        
        _zeta=np.rad2deg(-abs(_z)*np.sin(Omega+phi))*3600000 # mas/sec
        _eta=np.rad2deg(abs(_z)*np.cos(Omega+phi)*np.tan(np.radians(zeta)))*3600000-wz # mas/sec
        
        zeta_list.append(_zeta)
        eta_list.append(_eta)
        
    #plotEtaZeta(t_list, eta_list, zeta_list)
    
    return t_list, zeta_list, eta_list

def getSmearing(zeta_list, eta_list):
    dt=0.001 # 1 millisecond
    #AC_integration=sum(map(abs,zeta_list))*dt
    #AL_integration=sum(map(abs,[e+60000 for e in eta_list]))*dt
    AC_integration=sum(zeta_list)*dt
    AL_integration=sum([e+60000 for e in eta_list])*dt
    
    
    return AC_integration, AL_integration

    
if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 14})
    
    #buildFortranRoutine_init()
    w=Wrapper('CalWrapper.json')
    fa=FrequencyAnalysis(w)
    ip=ImgProc()
    
    o=w.getCollection()[0]
    print o.day_int+o.day_float

    start,Omega0,nu0=getInitParams(w)
    
    
    ###### original #######
    #Omega0=3.0823632010724396
    #nu0=5.8195930723884430
    
    ###### shifted #######
    #Omega0=2.6645082056929823
    #nu0=5.8180006687463219
    
    ###### after reboot ######
    #Omega0=-3.8220379251870327
    #nu0=2.8375593131104324
    
    
    
    print Omega0,nu0
    
    print
    start_date=o.timestamp
    #plotting 24 hours parameters
    #buildFortranRoutine()
    obmt_rev=obmt2rev(start_date)
    tot_hours=(w.observations[-1].timestamp-o.timestamp)*1.0/(1000000)
    dt=milliseconds2days(tot_hours*1.0/w.n_transits)
    lines=callFortranRoutine(obmt_rev,dt=dt,tmax=milliseconds2days(tot_hours),Omega0=Omega0,nu0=nu0)
    t_list, zeta_list, eta_list=getImageSpeed(lines, w.fov, o.ccdStrip, o.ccdRow)
    
    ################################# measuring shift ###################################
    hr_list=[(t1-t_list[0])*24 for t1 in t_list]
    hr2_list=[(o.timestamp-(w.getCollection(1)[0].timestamp))*1.0/(1000*1000*1000*60*60) for o in w.observations if o.ACrate>=0]
    zeta2_list=[o.ACrate for o in w.observations if o.ACrate>=0]
    
    
    wi, ha = figure.figaspect(6.)
    fig_sm, ax_sm = plt.subplots(figsize=(ha,wi))
    # Twin the x-axis twice to make independent y-axes.
    axes = [ax_sm, ax_sm.twinx()]
    # Make some space on the right side for the extra y-axis.
    fig_sm.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(hr_list,zeta_list, color='steelblue')
    axes[0].tick_params(axis='y', colors='steelblue')
    axes[0].tick_params(axis='y')
    axes[0].set_ylabel(r"Theoretical $\zeta$' [mas/s]")
    
    #time_list2=[x for x in time_list]
    axes[1].plot(hr2_list,zeta2_list,color='indianred')
    axes[1].tick_params(axis='y', colors='indianred')
    axes[1].tick_params(axis='y')
    axes[1].set_ylabel('Empirical AC rate')
    #axes[1].set_ylim(-5,5)
    #plt.xlim(0,24)
    plt.savefig('empirical_vs_theoretical.png', bbox_inches='tight')
    plt.close('all')
    print "fatto"
    wdt, h = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(h,wdt))
    
    ax.set_ylabel(r"$\zeta$' [mas/s]")
    ax.set_xlabel('time [OBMT days after J2000.0]')
    ax.set_ylim(-200,200)

    sm_list=[]
    rate_list=[]
    px_list=[]
    time_list=[]
    _z_list=[]
    ii=0
    txt_file_cooccurrence = open(fa.txt_cooccurrence,'w')
    txt_file_moments = open(fa.txt_moments,'w')
    for o in w.getCollection()[1:]:
        if ii%100==0 and o.ACrate>=0:
            hour=(o.timestamp-(w.getCollection(1)[0].timestamp))*1.0/(1000*1000*1000*60*60)
            ot=o.integration_time
            obmt=o.timestamp
            
            tmax=obmt-start_date
            
            tmax=milliseconds2days(tmax*1.0/1000000)
            #dt=milliseconds2days(1000)
            start_rev=obmt2rev(start_date)
            lines=callFortranRoutine(start_rev,dt=tmax,tmax=tmax,Omega0=Omega0,nu0=nu0)
            init = lines[-1].split()
            Omega0=init[3]
            nu0=init[5]
            #print o.id,Omega0,nu0
            
            obmt_rev=obmt2rev(obmt)
            dt=milliseconds2days(1)
            tmax=milliseconds2days(ot)
            lines=callFortranRoutine(obmt_rev,dt=dt,tmax=tmax,Omega0=Omega0,nu0=nu0)
    
            #lines=out.splitlines()
            _t, _z, _e=getImageSpeed(lines, w.fov, o.ccdStrip, o.ccdRow)
            
            AC_sm,AL_sm=getSmearing(_z, _e)
            if AC_sm>0:
                AC_px=(AC_sm+176.8)/176.8
            else:
                AC_px=(AC_sm-176.8)/176.8
            #AC_px=(AC_sm+645.249)/176.8
            if AL_sm>0:
                AL_px=(AL_sm+58.9)/58.9
            else:
                AL_px=(AL_sm-58.9)/58.9
            #AL_px=(AL_sm+222.5)/58.9
            aspect=(AL_px)/(AC_px)
            
            '''
            new_im=ip.img_resizeProfileNSL(o.window, AC_sm=AC_sm, AL_sm=AL_sm)
            res=ip.img_detectROI(new_im, thr_factor=0.5, zero_padding=False)
            roi=res[0]
            
            h,wi=np.shape(roi)
            aspect_roi=1.0*wi/h
            
            fa.writeROImoments(roi, o, txt_file_moments)
            fa.writeROIcooccurrence(roi, o, txt_file_cooccurrence)
            #print aspect_roi
            '''
            #print ('%s ---- %s ---- %s ---- %s ----- %s ---- %s ---- %s')%(str(o.id),str(round(aspect_roi,3)),str(round(AC_px,3)),str(round(AL_px,3)),str(round(o.ACrate,3)),str(round(AC_sm,3)),str(round(AL_sm,3)))
            print o.id
            
            if o.ACrate>=0:
                _z_list.append(abs(_z[0]))
                sm_list.append(AL_sm)
                px_list.append(AL_px)
                rate_list.append(o.ACrate)
                time_list.append(hour)
            
            
            start=min(_t)
            stop=max(_t)
            ax.plot(_t,_z,color='green', lw=1.3)
            ax.fill_between(_t,_z, facecolor='green', alpha=0.5)
            
            start_date=obmt
            '''
            fig2, axes = plt.subplots(nrows=1, ncols=2)
            axes[0].imshow(ip.img_interpolateImage(o.window, 1800, 1200))
            axes[0].set_title('Original (interpolated)')
            axes[1].imshow(new_im)
            axes[1].set_title('Shifted with NSL (not shrinked yet)')
            fig.suptitle('Image #%i' %o.id)
            plt.savefig('NSL/NSLresult'+str(o.id)+'.png', bbox_inches='tight')
            plt.close()
            '''
            #print o.id
        ii+=1
        
    txt_file_cooccurrence.close()
    txt_file_moments.close()
    
    '''
    wi, ha = figure.figaspect(6.)
    fig, ax = plt.subplots(figsize=(ha,wi))

    # Twin the x-axis twice to make independent y-axes.
    axes = [ax, ax.twinx()]
    
    # Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(time_list,zeta_list, color='steelblue', lw=0.7)
    axes[0].tick_params(axis='y', colors='steelblue')
    axes[0].set_ylabel('smearing')
    
    #time_list2=[x for x in time_list]
    axes[1].plot(time_list,rate_list,color='indianred')
    axes[1].tick_params(axis='y', colors='indianred')
    axes[1].set_ylabel('AC rate')
    
    plt.savefig('smearingVSrate.eps')
    plt.close('all')
    '''
    
    '''
    ################################### AC SMEARING PLOT ####################################
    wi, ha = figure.figaspect(6.)
    fig_sm, ax_sm = plt.subplots(figsize=(ha,wi))
    # Twin the x-axis twice to make independent y-axes.
    axes = [ax_sm, ax_sm.twinx()]
    # Make some space on the right side for the extra y-axis.
    fig_sm.subplots_adjust(right=0.75)
    axes[0].set_xlabel('hours')
    axes[0].plot(time_list,sm_list, color='white')
    axes[0].tick_params(axis='y')
    axes[0].set_ylabel('AL smearing [mas]')
    
    #time_list2=[x for x in time_list]
    axes[1].plot(time_list,px_list,color='indianred')
    axes[1].tick_params(axis='y')
    axes[1].set_ylabel('AL smearing [pix]')
    #axes[1].set_ylim(-5,5)
    #plt.xlim(0,24)
    plt.savefig('smearing_pix_mas.png', bbox_inches='tight')
    plt.close('all')
    '''
            
    
    '''
    xcorr=np.correlate(_z_list,rate_list)
    
    
    ab=np.argmax(signal.correlate(_z_list,rate_list))
    ba=np.argmax(signal.correlate(rate_list,_z_list))
    print ab
    print ba
    plt.plot(signal.correlate(_z_list,rate_list))
    plt.savefig('xcorr.eps')
    plt.close('all')
    
    af = fft(_z_list)
    bf = fft(rate_list)
    c = ifft(af * conj(bf))
    plt.plot(abs(c))
    plt.savefig('xcorr2.eps')
    plt.close('all')

    plt.plot(periodic_corr(rate_list, _z_list))
    print np.argmin(periodic_corr(rate_list, _z_list))
    print time_list[np.argmax(periodic_corr(rate_list, _z_list))]
    plt.savefig('xcorr3.eps')
    plt.close('all')
    '''
    
    plt.savefig('ACrateIntervals1.png', bbox_inches='tight')
    plt.close()
    
    
    