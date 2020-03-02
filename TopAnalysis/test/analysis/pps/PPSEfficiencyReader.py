import ROOT
import sys

class PPSEfficiencyReader:
    
    """ 
    takes care of reading the efficency measurements to memory and retrieving the final efficiency correction 
    see details in https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsStripsEfficiencies
    """

    def __init__(self, fList, year=2017):

        self.allEffs={}

        for fIn in fList.split(','):
            baseDir='Strips/%d'%year if 'MultiTrack' in fIn else 'Pixel/%d'%year
            fIn=ROOT.TFile.Open(ROOT.gSystem.ExpandPathName(fIn))
            for k in fIn.Get(baseDir).GetListOfKeys():
                for kk in fIn.Get(baseDir+'/'+k.GetName()).GetListOfKeys():
                    hname=kk.GetName()
                    if '2D'in hname : continue
                    self.allEffs[hname]=kk.ReadObj()
                    self.allEffs[hname].SetDirectory(0)
            fIn.Close()

        print '[PPSEfficiencyReader] retrieved %d histograms'%len(self.allEffs)


    def getPPSEfficiency(self,era,xangle,xi,rp,isMulti=True, applyMultiTrack=False):

        sector=45 if rp<100 else 56
        eff,effUnc=1.0,0.0
        if isMulti:

            if applyMultiTrack:
                multiTrack=self.allEffs['h%dmultitrackeff_%s_avg_RP%d'%(sector,era,rp)]
                ieff = multiTrack.GetBinContent(1)
                eff *= ieff
        
            raddam=self.allEffs['h%d_%s_%d_1D'%(sector,era,xangle)]
            raddamUnc=self.allEffs['h%derrors_%s_%d_1D'%(sector,era,xangle)]
            ibin=raddam.FindBin(xi)
            ieff    = raddam.GetBinContent(ibin)
            if ieff>0:
                eff    *= ieff
                effUnc += (raddamUnc.GetBinError(ibin)/ieff)**2            
                
        else:

            # FIXME
            # pixels are not  fully available yet so assume 100% eff
            eff,effUnc=1.0,0.0

            
        effUnc=eff*ROOT.TMath.Sqrt(effUnc)

        return eff,effUnc


    def getProjectedFinalState(self,
                               pos_protons,stripPosEff,stripPosEffUnc,
                               neg_protons,stripNegEff,stripNegEffUnc,
                               sighyp):

        """
        sighyp is a number between 0 and 16 where the bits represent
        0b - number of pixels in negative side
        1b - number of multi in negative side
        0b - number of pixels in positive side
        1b - number of multi in positive side
        """
           
        ppsWgt,ppsWgtUnc=1.0,0.0
        
        #check how many multi are required for the signal hypothesis
        nPixNegInSigHyp   = ((sighyp>>0) & 0x1)
        nMultiNegInSigHyp = ((sighyp>>1) & 0x1)
        nPixPosInSigHyp   = ((sighyp>>2) & 0x1)
        nMultiPosInSigHyp = ((sighyp>>3) & 0x1)

        #impossible cases! a multiRP needs a pixel
        if nMultiNegInSigHyp==1 and nPixNegInSigHyp==0:
            ppsWgt, ppsWgtUnc = 0., 0.
            pos_protons=[[],[],[]]
            neg_protons=[[],[],[]]
            return pos_protons,neg_protons,ppsWgt,ppsWgtUnc
        if nMultiPosInSigHyp==1 and nPixPosInSigHyp==0:
            ppsWgt, ppsWgtUnc = 0., 0.
            pos_protons=[[],[],[]]
            neg_protons=[[],[],[]]
            return pos_protons,neg_protons,ppsWgt,ppsWgtUnc

        #check how many multi and pixels are available
        nMultiPos = min(1,len(pos_protons[0]))
        nPixPos   = min(1,len(pos_protons[1]))
        nMultiNeg = min(1,len(neg_protons[0]))
        nPixNeg   = min(1,len(neg_protons[1]))

        #number of pixels must match!
        if nPixPosInSigHyp!=nPixPos:
            ppsWgt, ppsWgtUnc = 0., 0.
            pos_protons[1]=[]
            return pos_protons,neg_protons,ppsWgt,ppsWgtUnc
        if nPixNegInSigHyp!=nPixNeg :
            ppsWgt, ppsWgtUnc = 0., 0.
            neg_protons[1]=[]
            return pos_protons,neg_protons,ppsWgt,ppsWgtUnc

        #do the migrations due to inefficiency

        #cases where 1 proton is to be found on positive side
        if nMultiPosInSigHyp==1:

            #proton was already there, apply survival probability
            if nMultiPos==1 and stripPosEff!=0.: 
                ppsWgt        *= stripPosEff
                ppsWgtUnc     *= stripPosEffUnc/stripPosEff

            elif nMultiPos==0:
                ppsWgt, ppsWgtUnc = 0., 0.
                pos_protons[0]=[]
                pos_protons[2]=[]

        #cases where no proton is to be found on the positive side
        else:

            pos_protons[0] = []
            pos_protons[2] = []

            #in case one had been reconstructed downeight by inefficiency probability
            if nMultiPos==1 and stripPosEff<1: 
                ppsWgt        *= (1-stripPosEff)
                ppsWgtUnc     *= stripPosEffUnc/(1-stripPosEff)

        #cases where 1 proton is to be found on negative side
        if nMultiNegInSigHyp==1:

            #proton was already there, apply survival probability
            if nMultiNeg==1 and stripNegEff>0: 
                ppsWgt        *= stripNegEff
                ppsWgtUnc     *= stripNegEffUnc/stripNegEff
            elif nMultiNeg==0:
                ppsWgt, ppsWgtUnc = 0., 0.
                neg_protons[0]=[]
                neg_protons[2]=[]

        #cases where no proton is to be found on the negative side
        else :

            neg_protons[0]=[]
            neg_protons[2]=[]

            #in case one had been reconstructed downeight by inefficiency probability
            if nMultiNeg==1 and stripNegEff!=1.: 
                ppsWgt        *= (1-stripNegEff)
                ppsWgtUnc     *= stripNegEffUnc/(1-stripNegEff)

        #finalize weight uncertainty
        ppsWgtUnc= ppsWgt*ROOT.TMath.Sqrt(ppsWgtUnc)

        #return final result
        return pos_protons,neg_protons,ppsWgt,ppsWgt



def main():

    ppEffReader=PPSEfficiencyReader(fList='test/analysis/pps/PreliminaryEfficiencies_October92019_1D2DMultiTrack.root')

    #test for different conditions
    for era in ['2017B','2017C','2017D','2017E','2017F']:
        for xangle in [120,130,140,150]:
            for rp in [3,103]:
                eff=ppEffReader.getPPSEfficiency(era,xangle,0.035,rp)
                print '%6s %d %3d %3.3f +/- %3.3f'%(era,xangle,rp,eff[0],eff[1])


if __name__ == "__main__":
    sys.exit(main())
