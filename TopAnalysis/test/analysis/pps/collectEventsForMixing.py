import pickle
import os
import sys
from random import shuffle
from collections import defaultdict
from generateBinnedWorkspace import VALIDLHCXANGLES

#create mixing manks per era and crossing angle
for era in 'BCDEF':

    era='2017'+era
    for xangle in VALIDLHCXANGLES:

        rpData = defaultdict(list)

        baseDir=sys.argv[1]
        fList=[f for f in os.listdir(baseDir+'/mixing/Chunks') if era in f]

        toCheck=[]
        print 'Checking and merging',len(fList),'files for',(era,xangle)
        for f in fList:

            if 'ZeroBias' in f : continue
            if 'DoubleEG' in f : continue
            if 'Photon' in f : continue

            fullf=os.path.join(baseDir+'/mixing/Chunks',f)

            try:
                with open(fullf,'r') as cache:
                    a=pickle.load(cache)
                for key,evList in a.items():
                    kera,kangle,kevtype=key 
                    if kera   != era    : continue
                    if kangle != xangle : continue
                    rpData[key] += evList
            except Exception as e:
                toCheck.append(f)

        print '\t total number of events for mixing'
        for key in rpData:
            print '\t',key,len(rpData[key])

        mixbank='%s/mixing/mixbank_%s_%d.pck'%(baseDir,era,xangle)
        print '\t writing mixing bank @',mixbank
        with open('mixbank.pck','w') as cache:
            pickle.dump(rpData,cache)
        os.system('mv mixbank.pck %s'%mixbank)

if len(toCheck)>0:
    print '-'*50
    print 'Found these files with errors'
    print toCheck
    print '-'*50
