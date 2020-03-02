import ROOT
import os
import sys
import re
from TopLJets2015.TopAnalysis.Plot import *
from generateBinnedWorkspace import SIGNALXSECS

url=sys.argv[1]
tkns=re.findall(r'\d+',os.path.basename(url))
mass,xangle=tkns[0],tkns[1]
print mass,xangle
fIn=ROOT.TFile.Open(url)
nopu=fIn.Get('mmass_eenopu')
pu=fIn.Get('mmass_eerpinhpur')
data=fIn.Get('data')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

#
# before and after pileup mixing/selection
#
p=Plot('mmass_{0}_xangle{1}_sig'.format(mass,xangle),com='13 TeV')
p.frameMin=0
p.add(nopu,  title='no pileup',  color=ROOT.kGray, isData=False,  spImpose=False, isSyst=False)
p.add(pu,    title='pileup',     color=1,          isData=False, spImpose=True, isSyst=False)
finalExtraText='#scale[0.8]{{m={0} GeV #alpha={1}#murad}}'.format(mass,xangle)
p.show(outDir='./', lumi=37500,extraText=finalExtraText)

#
# in-fiducial vs out-fiducial components
#
fiducialCuts='gencsi1>0.03 & gencsi1<0.13 && gencsi2>0.03 && gencsi2<0.16'
data.Draw("mmiss>>hout(50,0,2500)","ppsEff*wgt*(!(%s) && mmiss>0 && mixType==1)"%fiducialCuts,'goff')
hout=ROOT.gDirectory.Get('hout')
hout.SetDirectory(0)
hout.Scale(SIGNALXSECS[int(xangle)]*37500.)
hout.GetXaxis().SetTitle('Missing mass [GeV]')
hout.GetYaxis().SetTitle('Events')

data.Draw("mmiss>>hin(50,0,2500)","ppsEff*wgt*(%s && mmiss>0 && mixType==1)"%fiducialCuts,"goff")
hin=ROOT.gDirectory.Get('hin')
hin.SetDirectory(0)
hin.Scale(SIGNALXSECS[int(xangle)]*37500.)
hin.SetLineWidth(2)
hin.GetXaxis().SetTitle('Missing mass [GeV]')
hin.GetYaxis().SetTitle('Events')



p=Plot('mmass_{0}_xangle{1}_sigacc'.format(mass,xangle),com='13 TeV')
p.add(hin,  title='in accept.',  color=ROOT.kRed,   isData=False, spImpose=False, isSyst=False)
p.add(hout, title='out accept.', color=ROOT.kBlack, isData=False, spImpose=True,  isSyst=False)
finalExtraText='#scale[0.8]{{m={0} GeV #alpha={1}#murad}}'.format(mass,xangle)
p.show(outDir='./', lumi=37500, extraText=finalExtraText)

