'''
Created on 18 mar 2018

@author: jimmijamma
'''
from os import listdir,makedirs
from os.path import isfile,join
import shutil as sh

DIRPATH='/Users/.../fotogaia'
HOURS=3

if __name__ == '__main__':
    onlyfiles = [f for f in listdir(DIRPATH) if isfile(join(DIRPATH,f))]
    intlist=[]
    for s in onlyfiles:
        fn=s.split('.')[0]
        if fn==None or fn=='':
            pass
        else:
            intlist.append(long(fn))
    
    n_bins = int((max(intlist)-min(intlist))/(HOURS*60*60*1000*1000*1000))+1
    bins = range(n_bins)
    for b in bins:
        makedirs(DIRPATH+'/'+str(b))
    for o in intlist:
        for b in bins:
            if o-min(intlist) < (b+1)*HOURS*60*60*1000*1000*1000:
                sh.copy2(DIRPATH+"/"+str(o)+".png",DIRPATH+'/'+str(b))
                break

        