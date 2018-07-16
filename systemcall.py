'''
Created on 15 lug 2018

@author: jimmijamma
'''
from subprocess import call,Popen,PIPE, STDOUT
import time


def buildFortranRoutine():
    call(['gfortran' ,'calcscanrate.for', '-o', 'calcscanrate.out'])
    time.sleep(2)
    

def callFortranRoutine(obmt_rev,dt,tmax):
    p=Popen(['./calcscanrate.out'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    input_string=str(obmt_rev)+','+str(dt)+','+str(tmax)
    sys_out=p.communicate(input=input_string)[0]
    return sys_out

def obmt2rev(obmt):
    one_rev=216e11
    obmt_rev=obmt*1.0/one_rev
    return obmt_rev

def milliseconds2days(ms):
    days=ms*1.0/86400000
    print days
    return days
    
if __name__ == '__main__':
    buildFortranRoutine()
    obmt=134270613714739200
    obmt_rev=obmt2rev(obmt)
    
    dt=milliseconds2days(1)
    tmax=milliseconds2days(4420)
    out=callFortranRoutine(obmt_rev,dt=dt,tmax=tmax)
    print out
    lines=out.splitlines()
    print len(lines)
