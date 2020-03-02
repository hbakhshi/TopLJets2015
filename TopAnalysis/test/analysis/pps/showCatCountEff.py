import ROOT
import sys

fIn=ROOT.TFile.Open(sys.argv[1])

for k in fIn.Get('catcount').GetListOfKeys():
    h=k.ReadObj()
    norm=h.GetBinContent(1)

    print h.GetTitle()
    totalNew=0
    for i in range(5):
        eff=h.GetBinContent(i+2)/norm
        print '%15s %3.4f +/- %3.4f'%(h.GetXaxis().GetBinLabel(i+2),eff,h.GetBinError(i+2)/norm)
        if i>0: totalNew+=eff
    print 'Total new: %3.4f'%totalNew
    print '-'*50


