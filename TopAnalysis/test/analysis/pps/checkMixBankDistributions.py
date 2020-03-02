import MixedEventSummary
import sys
import pickle
import ROOT

with open(sys.argv[1],'r') as f:
    rpData=pickle.load(f)

csi={}
for key in rpData:

    era,xangle,cat=key
    csi[cat]=[ROOT.TH1F('cat%d_csi%d'%(cat,i),';#xi;Protons',20,0,0.2) for i in range(3)]
    nevts=len(rpData[key])

    print key,'has',nevts,'events'
    for i in range(nevts):
        ev=rpData[key][i]
        protons=ev.pos_protons

        for i in range(3):            
            if len(protons[i])>0:
                csi[cat][i].Fill(max(protons[i]))
        if i>1000: break

fOut=ROOT.TFile.Open('rpdata.root','RECREATE')
for key in csi:
    for h in csi[key]:
        h.SetDirectory(fOut)
        h.Write()
fOut.Close()
