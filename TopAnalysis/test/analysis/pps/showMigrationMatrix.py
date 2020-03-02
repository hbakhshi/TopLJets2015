import ROOT
import sys
import os
import re
from generateBinnedWorkspace import SIGNALXSECS

preTS2=14586.4464/41529.3
postTS2=1-preTS2


tag=sys.argv[1]
inDir=sys.argv[2]
files=[os.path.join(inDir,f) for f in os.listdir(inDir) if tag in f] 
m=None
for f in files:

    xangle=int(re.search('xangle_(\d+)',f).group(1))

    inF=ROOT.TFile.Open(f)
    im=inF.Get('sighyp_wgt')
    if m is None:
        m=im.Clone('migmatrix')
        m.Reset('ICE')
        m.SetDirectory(0)

    sf  = preTS2 if 'preTS2' in f else postTS2
    sf *= SIGNALXSECS[xangle]

    m.Add(im,sf)
        
    inF.Close()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPaintTextFormat('4.3f');
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
c.SetGridx()
c.SetGridy()
m.Draw('coltext')
m.GetZaxis().SetTitle('Probability')
m.Scale(1./m.Integral())
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
c.Modified()
c.Update()
c.RedrawAxis()
c.SaveAs('%s_mig.png'%tag)
