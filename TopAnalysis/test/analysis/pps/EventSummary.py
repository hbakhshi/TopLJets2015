import ROOT
from array import array
import sys

class EventSummary:

    """ Event summary for final analysis """

    def __init__(self):

        """ define all variables needed to store and instantiate them as arrays of different types """

        self.vars={'i':['cat','isOffZ','xangle','run','lumi','era','sighyp','nch','nvtx','mixType','PFMultSumHF'],
                   'l':['event'],
                   'f':['wgt',
                        'l1pt','l1eta','l2pt','l2eta',
                        'bosonm','bosonpt','bosoneta','bosony','acopl','costhetacs',
                        'njets','mpf','zjb','zj2b','rho',
                        'PFHtSumHF','PFPzSumHF','rfc',
                        'gen_pzpp','gen_pzwgtUp','gen_pzwgtDown','gencsi1','gencsi2']
                   }

        for pfix in ['','syst']:
            self.vars['f'] += [pfix+'protonCat']
            addVars='{0}ppsEff,{0}ppsEffUnc,{0}csi1,{0}csi2,{0}mpp,{0}ypp,{0}pzpp,{0}mmiss,{0}ymmiss'.format(pfix)
            self.vars['f'] += addVars.split(',')

        self.vars['f'] += ['mmissvup','mmissvdn']

        #instantiate the arrays
        for t in self.vars:
            for v in self.vars[t]:
                setattr(self,v,array( t, [0] ))

    def reset(self):
        
        """reset values to 0 everywhere"""

        for t in self.vars:
            for v in self.vars[t]:
                getattr(self,v)[0]=0
        
        
    def attachToTree(self,tree):

        """ loop over types/variables and add the corresponding branches """

        for t in self.vars:
            for v in self.vars[t]:
                tree.Branch( v,   getattr(self,v), '%s/%s'%(v,t.upper()) )


def main():
    print 'Defines EventSummary class'

if __name__ == "__main__":
    sys.exit(main())
