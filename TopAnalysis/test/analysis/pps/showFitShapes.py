import sys
import ROOT
import re
from compareOptimResults import showShapes
from prepareOptimScanCards import OPTIMLIST


def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    lumi=37.5
    bosonName=sys.argv[3]
    if bosonName=='g' : 
        bosonName='#gamma'
        lumi=2.76

    optimPt=int(re.search('optim_(\d+)',sys.argv[1]).group(1))-1
    cuts=OPTIMLIST[optimPt][2].split(',')

    showShapes(resultsDir=sys.argv[1],
               name='shapes',
               title='pp%sX'%bosonName,
               mass=int(sys.argv[2]),
               boson=sys.argv[3],
               r95=None,
               sig=None,
               lumi=lumi,
               plotData=False, #True,
               showPseudoData=False,
               showAllBkgs=True,  #False)
               subCatTitles=cuts)

if __name__ == "__main__":
    main()
