import pickle
import ROOT
import random
import copy
from random import shuffle

class EventMixingTool:

    def __init__(self, mixedRP,validAngles):

        """ Reads the event mixing data from a pickle file and builds a list of crossing angle probabilities """

        self.mixedRP=mixedRP
        self.xangleRelFracs={}

        try:
            #build the list of probabilities for the crossing angles in each era
            print '[EventMixingTool] defining relative crossing angle fractions'
            gammaScaleFactors={'2017B': [1.0881694555282593, 1.0961500406265259, 1.0834789276123047, 0.8419778347015381], 
                               '2017C': [1.3276927471160889, 1.2163710594177246, 0.9180033206939697, 0.6171205639839172], 
                               '2017D': [1.306527853012085, 1.0550585985183716, 0.9354883432388306, 0.8344970941543579],                                
                               '2017E': [1.4099575281143188, 0.97742760181427, 0.7880721688270569, 0.6549458503723145], 
                               '2017F': [1.4421682357788086, 0.7824893593788147, 0.7729737162590027, 0.6496632695198059]}

            for key in self.mixedRP:
                mix_era,mix_xangle,mix_evcat=key
                if not mix_xangle in validAngles: continue
                n=len(self.mixedRP[key])
                xangleKey=(mix_era,mix_evcat)
                if not xangleKey in self.xangleRelFracs:
                    self.xangleRelFracs[xangleKey]=ROOT.TH1F('xanglefrac_%s_%d'%xangleKey,'',len(validAngles),0,len(validAngles))
                xbin=(mix_xangle-120)/10+1
                self.xangleRelFracs[xangleKey].SetBinContent(xbin,n)

            for xangleKey in self.xangleRelFracs.keys():
                self.xangleRelFracs[xangleKey].Scale(1./self.xangleRelFracs[xangleKey].Integral())                
                if xangleKey[1]!=13*13: continue
                newXangleKey=(xangleKey[0],22)
                self.xangleRelFracs[newXangleKey]=self.xangleRelFracs[xangleKey].Clone('xanglefrac_%s_%d'%newXangleKey)
                for xbin in range(4):
                    self.xangleRelFracs[newXangleKey].SetBinContent(xbin+1,
                                                                    gammaScaleFactors[xangleKey[0]][xbin]*self.xangleRelFracs[newXangleKey].GetBinContent(xbin+1))

        except Exception as e:
            print '-'*50
            print e
            print 'Note: this error is probably ok if one hasn\'t yet defined the mixing bank'
            print '-'*50


    def isIdle(self):
        
        """flag if it didn't initialize correctly"""

        return False if self.mixedRP else True


    def getNew(self,evEra,beamXangle,isData,validAngles,mixEvCategs,orderByDecreasingCsi=True):

        """get new list of mixed protons from different event categories
        also returns a list of variables which be used for pileup discrimination"""

        mixed_pos_protons={}
        mixed_neg_protons={}
        mixed_pudiscr={}

        try:
            for mixEvCat in mixEvCategs:
                mixed_pos_protons[mixEvCat]=[[],[],[]]
                mixed_neg_protons[mixEvCat]=[[],[],[]]
                mixed_pudiscr[mixEvCat]=[]

                if isData and not beamXangle in validAngles : continue
                
                mixedEvKey                  = (evEra,beamXangle,mixEvCat)
                mixedEv                     = random.choice( self.mixedRP[mixedEvKey] )                
                mixed_pudiscr[mixEvCat]     = copy.deepcopy(mixedEv.puDiscr)
                mixed_pos_protons[mixEvCat] = copy.deepcopy(mixedEv.pos_protons)               
                mixed_neg_protons[mixEvCat] = copy.deepcopy(mixedEv.neg_protons)

                if orderByDecreasingCsi:
                    mixed_pos_protons[mixEvCat][1].sort(reverse = True) 
                    mixed_neg_protons[mixEvCat][1].sort(reverse = True) 

        except Exception as e:
            print e  
            print evEra,beamXangle,mixEvCat,'->',mixedEvKey
            pass

        return mixed_pos_protons, mixed_neg_protons, mixed_pudiscr


    def mergeWithMixedEvent(self,
                            pos_protons, mixed_pos_protons,
                            neg_protons, mixed_neg_protons,
                            orderByDecreasingCsi = True):
            
        """merges tracks from two different events (useful for pure signal embedding) """

        merged_pos_protons={}
        merged_neg_protons={}

        #iterate over sources of mixed events
        for mixEvCat in mixed_pos_protons:

            merged_pos_protons[mixEvCat]=[[],[],[]]
            merged_neg_protons[mixEvCat]=[[],[],[]]

            #iterate over track reconstruction methods
            #for multiRP and strips if >1 track, reset (only 1 allowed per event in 2017)
            for i in range(3):

                merged_pos_protons[mixEvCat][i] += mixed_pos_protons[mixEvCat][i]
                merged_pos_protons[mixEvCat][i] += pos_protons[i]

                merged_neg_protons[mixEvCat][i] += mixed_neg_protons[mixEvCat][i]
                merged_neg_protons[mixEvCat][i] += neg_protons[i]

                if orderByDecreasingCsi:
                    merged_pos_protons[mixEvCat][i].sort(reverse = True) 
                    merged_neg_protons[mixEvCat][i].sort(reverse = True) 
                else:
                    shuffle(merged_pos_protons[mixEvCat][i])                
                    shuffle(merged_neg_protons[mixEvCat][i])

                if i==1: continue

                if len(merged_pos_protons[mixEvCat][i])>1:
                    merged_pos_protons[mixEvCat][i]=[]

                    #make sure multiRP have been killed by > 1 strip
                    if i==2: merged_pos_protons[mixEvCat][0]=[]

                if len(merged_neg_protons[mixEvCat][i])>1:
                    merged_neg_protons[mixEvCat][i]=[]

                    #make sure multiRP have been killed by > 1 strip
                    if i==2: merged_neg_protons[mixEvCat][0]=[]

        return merged_pos_protons, merged_neg_protons


    def getRandomLHCCrossingAngle(self,evEra,evCat):

        """samples the LHC crossing angle according to the relative fractions"""
        
        xbin=-1
        if self.xangleRelFracs:
            xangleKey=(evEra,evCat)
            xbin=ROOT.TMath.FloorNint(self.xangleRelFracs[xangleKey].GetRandom())

        return xbin
