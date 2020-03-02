import ROOT

from TopLJets2015.TopAnalysis.Plot import *
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sig='/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/Chunks/MC13TeV_ZX960_fullsim_0.root'
sigF=ROOT.TFile.Open(sig)
sigT=sigF.Get('tree')

dyT=ROOT.TChain('tree')
for i in range(5):
    dyT.AddFile('/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/Chunks/MC13TeV_2017_DY50toInf_fxfx_%d.root'%i)


finalSel='bosonpt>50 && isZ && abs(l1id*l2id)==169 && (hasMMTrigger || hasMTrigger)'
for t,v,c in [("p_{T}(ll) [GeV]",                 "bosonpt >> {0}(20,0,250)",   'isZ'),
              ('m(ll) [GeV]',                     "mboson >> {0}(20,81,101)",    finalSel),
              ('Jet multiplicity',                "nj >> {0}(4,0,4)",            finalSel),              
              ('Missing transverse energy [GeV]', "met_pt >> {0}(20,0,200)",     finalSel),
              ('Charged particle multiplicity',   "nchPV >> {0}(20,0,100)",       finalSel),
              ('p_{T}(charged)',                  'sumPVChPt >> {0}(20,0,250)',  finalSel),
              ('N(muons)',                        'nrawmu-2 >> {0}(10,0,10)',      finalSel),
              ]:
    sigT.Draw(v.format('sig'),'evwgt*(%s)'%c)
    sigH=ROOT.gDirectory.Get('sig')
    sigH.Scale(1./sigH.Integral())
    sigH.GetXaxis().SetTitle(t)
    sigH.GetYaxis().SetTitle('PDF')

    dyT.Draw(v.format('dy'),'evwgt*(%s)'%c)
    dyH=ROOT.gDirectory.Get('dy')
    dyH.Scale(1./dyH.Integral())
    dyH.GetXaxis().SetTitle(t)
    dyH.GetYaxis().SetTitle('PDF')

    p=Plot(v.split(' >>')[0],com='13 TeV')
    p.frameMin=0
    p.add(dyH,   title='DY',                                 color=ROOT.kMagenta+1, isData=False,  spImpose=True, isSyst=False)
    p.add(sigH,  title='#splitline{m_{X}=950 GeV}{#scale[0.5]{120#murad post-TS2}}',  color=ROOT.kBlack,     isData=False,  spImpose=True, isSyst=False)
    p.show(outDir='./', lumi=37500,)
    sigH.Delete()
    dyH.Delete()
    
