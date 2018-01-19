import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle
import numpy
from subprocess import Popen, PIPE, STDOUT

from TopLJets2015.TopAnalysis.Plot import *
from TopLJets2015.TopAnalysis.dataCardTools import *

def noPartIn(a,b) :
        return 0 == len([c for c in b if (a in c or c in a)])

def buildRateUncertainty(varDn,varUp):
    """ returns a uniformized string for the datacard for the treatment of up/down rate uncertainties """

    #if this is below 0.1% neglect it
    if abs(varUp-1)<0.001 and abs(varDn-1)<0.001: return None

    #distinguish same sided from double sided variations
    toReturn=None
    if varUp>1. and varDn>1.:
        satUnc=min(max(varUp,varDn),2.0)
        toReturn='%3.3f'%satUnc
    elif varUp<1. and varDn<1.:
        satUnc=max(min(varUp,varDn),1/2.0)
        toReturn='%3.3f'%satUnc
    else:
        toReturn='%3.3f/%3.3f'%(varDn,varUp)

    #all done here
    return toReturn

def getDistsFromDirIn(url,indir,applyFilter=''):
    """customize getting the distributions for hypothesis testing"""
    fIn=ROOT.TFile.Open(url)
    obs,exp=getDistsFrom(fIn.Get(indir),applyFilter)
    fIn.Close()
    return obs,exp

def getRowFromTH2(tempHist2D,columnName) :
    tempBinNum = tempHist2D.GetYaxis().FindBin(columnName);
    new1DProj  = tempHist2D.ProjectionX(tempHist2D.GetName()+"_"+columnName,
            tempBinNum,
            tempBinNum).Clone(tempHist2D.GetName())
    new1DProj.SetDirectory(0)

    return new1DProj


def getDistsForHypoTest(cat,rawSignalList,opt,outDir="",systName="",systIsGen=False):
    """
    readout distributions from ROOT file and prepare them to be used in the datacard
    """
    # do we want systematics?
    systExt = ""
    if systIsGen : systExt = "_gen"
    elif systName != "" : systExt = "_exp"

    # get main dists
    obs,exp=getDistsFromDirIn(opt.input,'%s_%s_w100%s'%(cat,opt.dist,systExt))

    #run through exp and convert into TH1s
    if systName != "" :
        for tHist2D in exp :
            exp[tHist2D] = getRowFromTH2(exp[tHist2D],systName)

    #signal hypothesis
    expMainHypo=exp.copy()
    if opt.mainHypo!=100.0:
        _,expMainHypo=getDistsFromDirIn(opt.input,'%s_%s_w%.0f'%(cat,opt.dist,opt.mainHypo))
    expAltHypo=None
    if len(opt.altHypoFromSim)>0 :
        _,expAltHypo=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f%s'%(cat,opt.dist,opt.altHypo,systExt))
        expAltHypo = {k.replace(opt.altHypoFromSim,""): v for k, v in expAltHypo.items() if opt.altHypoFromSim in v.GetName()}
        if systName != "" :
            for tHist2D in expAltHypo :
                expAltHypo[tHist2D] = getRowFromTH2(expAltHypo[tHist2D],systName)
    else:
        _,expAltHypo=getDistsFromDirIn(opt.input,'%s_%s_w%.0f%s'%(cat,opt.dist,opt.altHypo,systExt))
        if systName != "" :
            for tHist2D in expAltHypo :
                expAltHypo[tHist2D] = getRowFromTH2(expAltHypo[tHist2D],systName)

    #replace DY shape from alternative sample
    try:
        if opt.replaceDYshape:
            _,altDY=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f'%(cat,opt.dist,opt.mainHypo),'DY')
            nbins=exp['DY'].GetNbinsX()
            sf=exp['DY'].Integral(0,nbins+1)/altDY['DY'].Integral(0,nbins+1)
            altDY['DY'].Scale(sf)

            #save a plot with the closure test
            if outDir!="" and opt.doValidation:
                try:
                    plot=Plot('DY_%s_%s'%(cat,opt.dist))
                    plot.savelog=False
                    plot.doChi2=True
                    plot.wideCanvas=False
                    plot.plotformats=['pdf','png']
                    plot.add(exp['DY'],  "MG5_aMC@NLO FxFx (NLO)", 1, False, False)
                    plot.add(altDY['DY'],"Madgraph MLM (LO)",      2, False, False)
                    plot.finalize()
                    plot.show(outDir=outDir,lumi=12900,noStack=True)
                except:
                    pass

            #all done here
            exp['DY'].Delete()
            exp['DY']=altDY['DY']
    except:
        pass

    # add signal hypothesis to expectations
    # we normalize the expectations to the standard expectation
    # as we are looking for shape distortions due to the width
    for proc in rawSignalList:
        try:
            newProc=('%sw%.0f'%(proc,opt.mainHypo)).replace('.','p')

            #main hypothesis
            exp[newProc]=expMainHypo[proc].Clone(newProc)
            nbins=exp[newProc].GetNbinsX()
            sf=exp[proc].Integral(0,nbins+1)/exp[newProc].Integral(0,nbins+1)
            exp[newProc].Scale(sf)
            exp[newProc].SetDirectory(0)

            #alternative hypothesis
            newProc=('%sw%.0f'%(proc,opt.altHypo)).replace('.','p')
            if opt.mainHypo==opt.altHypo: newProc+='a'
            exp[newProc]=expAltHypo[proc].Clone(newProc)
            sf=exp[proc].Integral(0,nbins+1)/exp[newProc].Integral(0,nbins+1)
            exp[newProc].Scale(sf)
            exp[newProc].SetDirectory(0)
        except:
            pass

    #delete the nominal expectations
    for proc in rawSignalList:
        try:
            exp[proc].Delete()
            del exp[proc]
        except:
            pass

    return obs,exp


"""
prepare the steering script for combine
"""
def doCombineScript(opt,args,outDir,dataCardList):

    altHypoTag=('w%.0f'%opt.altHypo).replace('.','p')
    if opt.altHypo==opt.mainHypo : altHypoTag+='a'

    scriptname='%s/steerHypoTest.sh'%outDir
    script=open(scriptname,'w')
    print 'Starting script',scriptname
    script.write('#\n')
    script.write('# Generated by %s with git hash %s for standard (alternative) hypothesis %.0f (%.0f)\n' % (getpass.getuser(),
                                                                                                               commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                                                                                                               opt.mainHypo,
                                                                                                               opt.altHypo) )
    script.write('### environment setup\n')
    script.write('COMBINE=%s\n'%opt.combine)
    script.write('SCRIPTDIR=`dirname ${0}`\n')
    script.write('cd ${COMBINE}\n')
    script.write('eval `scramv1 r -sh`\n')
    script.write('cd ${SCRIPTDIR}\n')
    script.write('\n')

    script.write('### combine datacard and start workspace\n')
    script.write('combineCards.py %s > datacard.dat\n'%dataCardList)
    script.write('\n')

    script.write('### convert to workspace\n')
    script.write('text2workspace.py datacard.dat -P HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest -m 172.5 --PO verbose --PO altSignal=%s --PO muFloating -o workspace.root\n'%altHypoTag)
    script.write('\n')

    if opt.doValidation:
        script.write('### dump systematics\n')
        script.write('python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py datacard.dat --all -f html > systs_summary.html;\n')
        script.write('\n')

        script.write('### likelihood scans and fits\n')
        commonOpts="-m 172.5 --setPhysicsModelParameters x=${x},r=1 --setPhysicsModelParameterRanges r=0.8,1.2 --saveWorkspace --robustFit 1 --minimizerAlgoForMinos Minuit2,Migrad --minimizerAlgo Minuit2"
        #script.write('for x in 0 1; do\n')
        #script.write('   combine workspace.root -M MultiDimFit -P x --floatOtherPOI=1  --algo=grid --points=50 -t -1 --expectSignal=1 -n x_scan_${x}_exp %s;\n'%commonOpts)
        #script.write('   combine workspace.root -M MaxLikelihoodFit --redefineSignalPOIs x -t -1 --expectSignal=1 -n x_fit_${x}_exp --saveWithUncertainties %s;\n'%commonOpts)
        #script.write('done\n')
        script.write('combine workspace.root -M MultiDimFit -P x --floatOtherPOI=1 --algo=grid --points=50 -n x_scan_obs %s\n'%commonOpts)
        script.write('combine workspace.root -M MultiDimFit -P r --floatOtherPOI=1 --algo=grid --points=50 -n r_scan_obs %s\n'%commonOpts)
        script.write('combine workspace.root -M MaxLikelihoodFit --saveWithUncertainties --redefineSignalPOIs x -n x_fit_obs %s;\n'%commonOpts)
        script.write('combine workspace.root -M MaxLikelihoodFit --saveWithUncertainties --redefineSignalPOIs r -n r_fit_obs %s;\n'%commonOpts)
        script.write('\n')
        #these lines are commented out but can be useful for further local debugging
        commonOpts="--minimizerAlgo Minuit2 --minimizerStrategy 2 --algo=saturated -m 172.5  --redefineSignalPOIs r --setPhysicsModelParameterRanges r=0.8,1.2"
        script.write('#combine -M GoodnessOfFit workspace.root %s\n'%commonOpts)
        script.write('#combine -M GoodnessOfFit workspace.root %s -t 100 --fixedSignalStrength=1\n'%commonOpts)
        script.write('#mv higgsCombineTest.GoodnessOfFit.mH172.5.root        gof_r_obs.root\n')
        script.write('#mv higgsCombineTest.GoodnessOfFit.mH172.5.123456.root gof_r_exp.root\n')
        script.write('\n')
        commonOpts="--minimizerAlgo Minuit2 --minimizerAlgoForMinos Minuit2,Migrad --robustFit 1 -m 172.5 --redefineSignalPOIs r --setPhysicsModelParameterRanges r=0.8,1.2"
        script.write('#combineTool.py -M Impacts -d workspace.root %s --doInitialFit\n'%commonOpts)
        script.write('#combineTool.py -M Impacts -d workspace.root %s --doFits\n'%commonOpts)
        script.write('#combineTool.py -M Impacts -d workspace.root %s -o impacts.json\n'%commonOpts)
        script.write('#plotImpacts.py -i impacts.json -o impacts\n')
        script.write('\n')

    script.write('### SCAN \n')
    script.write('\n')
    for extra,extraName in [('-S 0','_stat'),('','')]:
        commonOpts="-m 172.5 %s -M HybridNew --testStat=TEV --onlyTestStat --saveToys --saveHybridResult --minimizerAlgo Minuit2"%extra
        if opt.frzString != "" :
            commonOpts += " --freezeNuisances %s"%opt.frzString
        script.write("combine %s --singlePoint 0  workspace.root -n scan0n\n"%commonOpts)
        script.write("mv higgsCombinescan0n.HybridNew.mH172.5.123456.root testStat_scan0n%s.root\n"%extraName)
        script.write("combine %s --singlePoint 1  workspace.root -n scan1n\n"%commonOpts)
        script.write("mv higgsCombinescan1n.HybridNew.mH172.5.123456.root testStat_scan1n%s.root\n"%extraName)


    #script.write('### CLs\n')
    # do not write CLs -- python can't launch scripts with forking
    #script.write('combine workspace.root -M HybridNew --seed 8192 --saveHybridResult -m 172.5 --saveWorkspace --saveToys --testStat=TEV --singlePoint 1 -T %d -i 2 --fork 6 --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 --expectedFromGrid 0.5 -n cls_prefit_exp;\n'%opt.nToys)
    #script.write('combine workspace.root -M HybridNew --seed 8192 --saveHybridResult -m 172.5 --saveWorkspace --saveToys --testStat=TEV --singlePoint 1 -T %d -i 2 --fork 6 --clsAcc 0 --fullBToys  --frequentist --expectedFromGrid 0.5 -n cls_postfit_exp;\n'%opt.nToys)
    #script.write('combine workspace.root -M HybridNew --seed 8192 --saveHybridResult -m 172.5 --saveWorkspace --saveToys --testStat=TEV --singlePoint 1 -T %d -i 2 --fork 6 --clsAcc 0 --fullBToys  --frequentist -n cls_postfit_obs;\n'%opt.nToys)
    script.write('\n')
    script.close()

    return scriptname

"""
instantiates one datacard per category
"""
def doDataCards(opt,args):

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')
    ttScenarioList=['tbart']
    mainSignalList,altSignalList=[],[]
    if 'tbart' in rawSignalList:
        ttScenarioList = [('tbartw%.0f'%h).replace('.','p') for h in [opt.mainHypo,opt.altHypo]]
        if opt.mainHypo==opt.altHypo: ttScenarioList[1]+='a'
        mainSignalList += [ttScenarioList[0]]
        altSignalList  += [ttScenarioList[1]]
    tWScenarioList=['Singletop']
    if 'Singletop' in rawSignalList:
        tWScenarioList = [('Singletopw%.0f'%h).replace('.','p') for h in [opt.mainHypo,opt.altHypo]]
        if opt.mainHypo==opt.altHypo: tWScenarioList[1]+='a'
        mainSignalList += [tWScenarioList[0]]
        altSignalList  += [tWScenarioList[1]]

    if opt.rmvNuisances != "" :
        rmvNuisances = True
        nuisanceRMV = opt.rmvNuisances.split(',')
    else :
        rmvNuisances = False

    if opt.frzNuisances != "" :
        frzNuisances = True
        nuisanceFRZ = opt.frzNuisances.split(',')
    else :
        frzNuisances = False


    #define RATE systematics : syst,val,pdf,whiteList,blackList  (val can be a list of values [-var,+var])
    rateSysts=[
          ('lumi_13TeV',       1.025,    'lnN',    [],                  ['DY','W']),
          ('DYnorm_*CH*',      1.30,     'lnN',    ['DY'],              []),
          ('Wnorm_th',         1.50,     'lnN',    ['W'],               []),
          ('tWnorm_th',        1.15,     'lnN',    tWScenarioList,      []),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson'],      []),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV'],          []),
    ]

    #define the SHAPE systematics from weighting, varying object scales, efficiencies, etc.
    # syst,weightList,whiteList,blackList,shapeTreatement=0 (none), 1 (shape only), 2 (factorizeRate),nsigma
    # a - in front of the process name in the black list will exclude rate uncertainties
    weightingSysts=[
        ('ees',            ['ees'],                                    [],             ['DY','-W'], 2, 1.0),
        ('mes',            ['mes'],                                    [],             ['DY','-W'], 2, 1.0),
        ('jer',            ['jer'],                                    [],             ['DY','-W'], 2, 1.0),
        ('trig_*CH*',      ['trig'],                                   [],             ['DY','-W'], 2, 1.0),
        ('sel_E',          ['esel'],                                   [],             ['DY','-W'], 2, 1.0),
        ('sel_M',          ['msel'],                                   [],             ['DY','-W'], 2, 1.0),
        ('ltag',           ['ltag'],                                   [],             ['DY','-W'], 2, 1.0),
        ('btag',           ['btag'],                                   [],             ['DY','-W'], 2, 1.0),
        ('bfrag',          ['bfrag'],                                  [],             ['DY','-W'], 2, 1.0),
        ('semilep',        ['semilep'],                                [],             ['DY','-W'], 2, 1.0),
        ('pu',             ['pu'],                                     [],             ['DY','-W'], 1, 1.0),
        ('tttoppt',        ['toppt'],                                  ttScenarioList, [],          2, 1.0),
        ('ttMEqcdscale',   ['gen%d'%ig for ig in[3,5,6,4,8,10] ],      ttScenarioList, [],          1, 1.0),
        ('ttPDF',          ['gen%d'%(11+ig) for ig in xrange(0,100) ], ttScenarioList, [],          0, 1.0)
        ]
    for ig in xrange(0,29) :
        weightingSysts += [('jes%s'%ig,            ['jes%d'%ig],       [],             ['DY'], 2, 1.0)]

    #define the SHAPE systematics from dedicated samples : syst,{procs,samples}, shapeTreatment (see above) nsigma
    fileShapeSysts = [
        ('mtop',           {'tbart':['t#bar{t} m=171.5',  't#bar{t} m=173.5']}       , 1, 1./2.),
        ('st_wid',         {'Singletop':['Single top m=169.5', 'Single top m=175.5']}, 1, 1./6.),
        ('UE',             {'tbart':['t#bar{t} UEdn',     't#bar{t} UEup']}          , 2, 1.0 ),
        ('CR',             {'tbart':['t#bar{t} QCDbased', 't#bar{t} gluon move']}    , 2, 1.0 ),
        ('hdamp',          {'tbart':['t#bar{t} hdamp dn', 't#bar{t} hdamp up']}      , 2, 1.0 ),
        ('ISR_tt',         {'tbart':['t#bar{t} isr dn',   't#bar{t} isr up']}        , 2, 1.0 ),
        ('FSR_tt',         {'tbart':['t#bar{t} fsr dn',   't#bar{t} fsr up']}        , 2, 1.0 ),
        ('ISR_st',         {'Singletop':['Single top isr dn', 'Single top isr up']}  , 2, 1.0 ),
        ('FSR_st',         {'Singletop':['Single top fsr dn', 'Single top fsr up']}  , 2, 1.0 ),
        ('tWttInterf',     {'Singletop':   ['Single top DS']}                        , 2, 1.0 ),
        ('tWMEScale',      {'Singletop':   ['Single top me dn', 'Single top me up']} , 2, 1.0 ),
        ]

    if rmvNuisances:
        rateSysts      = [a for a in rateSysts      if noPartIn(a[0],nuisanceRMV)]
        weightingSysts = [a for a in weightingSysts if noPartIn(a[0],nuisanceRMV)]
        fileShapeSysts = [a for a in fileShapeSysts if noPartIn(a[0],nuisanceRMV)]
        print "\n"
        print rateSysts
        print "\n"
        print weightingSysts
        print "\n"
        print fileShapeSysts

    # really convoluted, but this was the best way, I promise
    if frzNuisances and "all" not in nuisanceFRZ:

        # collect all correct systematic names, including %sRate
        frzString=",".join([",".join([a[0] for a in rateSysts      if not noPartIn(a[0],nuisanceFRZ)]),
                            ",".join([a[0] for a in weightingSysts if not noPartIn(a[0],nuisanceFRZ)]),
                            ",".join([a[0] for a in fileShapeSysts if not noPartIn(a[0],nuisanceFRZ)])])
        frzString+=","
        frzString+=",".join([",".join([a[0]+"Rate" for a in weightingSysts if not noPartIn(a[0],nuisanceFRZ) and a[4]==2]),
                             ",".join([a[0]+"Rate" for a in fileShapeSysts if not noPartIn(a[0],nuisanceFRZ) and a[2]==2])])
        frzString=frzString.replace(',,,',',')
        frzString=frzString.replace(',,',',')
        frzString=frzString[1:] if frzString[0] == "," else frzString
        frzString=frzString[:-1] if frzString[-1] == "," else frzString

        # add in more if splitting by channel
        if "*CH*" in frzString :
            catList = opt.cat.split(',')
            frzString = frzString.split(',')

            for elem in frzString :
                if "*CH*" in elem :
                    newelem=[]
                    if len([x for x in catList if "EE" in x]):
                        newelem += [elem.replace('*CH*','EE')]
                    if len([x for x in catList if "EM" in x]):
                        newelem += [elem.replace('*CH*','EM')]
                    if len([x for x in catList if "MM" in x]):
                        newelem += [elem.replace('*CH*','MM')]

                    newelem   = ','.join(newelem)
                    frzString.remove(elem)
                    frzString += [newelem]

            frzString = ','.join(frzString)

        opt.frzString=frzString
        print "\n"
        print frzString
        print "\n"
    elif "all" in nuisanceFRZ :
        frzString = "all"

        opt.frzString=frzString
        print "\n"
        print frzString
        print "\n"


    # prepare output directory
    outDir='%s/hypotest_%.0fvs%.0f%s'%(opt.output, opt.mainHypo,opt.altHypo,'sim'+opt.altHypoFromSim if len(opt.altHypoFromSim)!=0 else '')
    if opt.pseudoData==-1 : outDir += '_data'
    else:
        outDir += '_%.0f'%opt.pseudoData
        if len(opt.pseudoDataFromSim)!=0   : outDir+='sim_'
        elif len(opt.pseudoDataFromWgt)!=0 : outDir+='wgt_'
        outDir += 'pseudodata'
    os.system('mkdir -p %s'%outDir)
    os.system('rm -rf %s/*'%outDir)

    # prepare output ROOT file
    outFile='%s/shapes.root'%outDir
    fOut=ROOT.TFile.Open(outFile,'RECREATE')
    fOut.Close()

    # parse the categories to consider
    dataCardList=''
    for cat in opt.cat.split(','):
        lfs='EE'
        if 'EM' in cat : lfs='EM'
        if 'MM' in cat : lfs='MM'

        #data and nominal shapes
        obs,exp=getDistsForHypoTest(cat,rawSignalList,opt,outDir)

        #recreate data if requested
        if opt.pseudoData!=-1:
            pseudoSignal=None
            print '\t pseudo-data is being generated',
            if len(opt.pseudoDataFromSim) and opt.systInput:
                print 'injecting signal from',opt.pseudoDataFromSim
                pseudoDataFromSim=opt.pseudoDataFromSim.replace('_',' ')
                _,pseudoSignalRaw=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f'%(cat,opt.dist,opt.mainHypo),pseudoDataFromSim)
                pseudoSignal={}
                pseudoSignal['tbart']=pseudoSignalRaw.popitem()[1]
            elif len(opt.pseudoDataFromWgt):
                print 'injecting signal from',opt.pseudoDataFromWgt
                _,pseudoSignal=getDistsFromDirIn(opt.input,'%s%s_%s_w%.0f'%(opt.pseudoDataFromWgt,cat,opt.dist,opt.mainHypo),'t#bar{t}')
                print pseudoSignal,'%s%s_%s_w%.0f'%(opt.pseudoDataFromWgt,cat,opt.dist,opt.mainHypo)
            else:
                print 'injecting signal from weighted',opt.pseudoData
                _,pseudoSignal=getDistsFromDirIn(opt.input,'%s_%s_w%.0f'%(cat,opt.dist,opt.pseudoData))
            obs.Reset('ICE')

            #build pseudo-expectations
            pseudoSignalAccept=[]
            for proc in pseudoSignal:
                accept=False
                for sig in rawSignalList:
                    if sig==proc: accept=True
                if not accept : continue
                print "\t\t Including:", proc, pseudoSignal[proc].GetName()

                newProc=('%sw%.0f'%(proc,opt.mainHypo)).replace('.','p')
                pseudoSignalAccept.append(newProc)
                sf=exp[newProc].Integral()/pseudoSignal[proc].Integral()
                pseudoSignal[proc].Scale(sf)
                if opt.rndmPseudoSF :
                    from random import uniform
                    pseudoSignal[proc].Scale(uniform(0.99,1.01))
                obs.Add( pseudoSignal[proc] )

            if len(opt.pseudoDataFromWgt) : pseudoSignalAccept+=altSignalList

            for proc in exp:
                if "%.0f"%opt.altHypo in proc : continue
                if not proc in pseudoSignalAccept:
                    print "\t\t Including:", proc, exp[proc].GetName()
                    obs.Add( exp[proc] )
            print pseudoSignalAccept
            for xbin in xrange(0,obs.GetNbinsX()+2): obs.SetBinContent(xbin,int(obs.GetBinContent(xbin)))

        #start the datacard header
        datacardname='%s/datacard_%s.dat'%(outDir,cat)
        dataCardList+='%s=%s '%(cat,os.path.basename(datacardname))
        datacard=open(datacardname,'w')
        print 'Starting datacard',datacardname
        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis category %s\n' % (getpass.getuser(),
                                                                                          commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                                                                                          cat) )

        datacard.write('#\n')
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')
        datacard.write('-'*50+'\n')
        datacard.write('shapes *        * shapes.root %s_%s/$PROCESS %s_%s_$SYSTEMATIC/$PROCESS\n'%(cat,opt.dist,cat,opt.dist))

        #observation
        datacard.write('-'*50+'\n')
        datacard.write('bin 1\n')
        datacard.write('observation %3.1f\n' % obs.Integral())

        #nominal expectations
        print '\t nominal expectations',len(exp)-1
        datacard.write('-'*50+'\n')
        datacard.write('\t\t\t %16s'%'bin')
        for i in xrange(0,len(exp)): datacard.write('%15s'%'1')
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'process')
        for sig in mainSignalList: datacard.write('%15s'%sig)
        for sig in altSignalList:  datacard.write('%15s'%sig)
        for proc in exp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'process')
        procCtr=-len(mainSignalList)-len(altSignalList)+1
        for sig in mainSignalList:
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        for sig in altSignalList:
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        for proc in exp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'rate')
        for sig in mainSignalList: datacard.write('%15s'%('%3.2f'%(exp[sig].Integral())))
        for sig in altSignalList:
            #if 'Singletop' in sig : sig = 'Singletopw100'
            datacard.write('%15s'%('%3.2f'%(exp[sig].Integral())))
        for proc in exp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%('%3.2f'%(exp[proc].Integral())))
        datacard.write('\n')
        datacard.write('-'*50+'\n')

        #save to nominal to shapes file
        nomShapes=exp.copy()
        nomShapes['data_obs']=obs
        #for h in exp: nomShapes[h]=exp[h].Clone( h+"_final" )
        #nomShapes['data_obs']=obs.Clone('data_obs_final')
        saveToShapesFile(outFile,nomShapes,('%s_%s'%(cat,opt.dist)),opt.rebin)

        #MC stats systematics for bins with large stat uncertainty
        if opt.addBinByBin>0:
            for proc in exp:
                finalNomShape=exp[proc].Clone('tmp')
                if opt.rebin : finalNomShape.Rebin(opt.rebin)

                for xbin in xrange(1,finalNomShape.GetXaxis().GetNbins()+1):
                    val,unc=finalNomShape.GetBinContent(xbin),finalNomShape.GetBinError(xbin)
                    if val==0 : continue
                    if ROOT.TMath.Abs(unc/val)<opt.addBinByBin: continue

                    binShapes={}
                    systVar='%sbin%d%s'%(proc,xbin,cat)

                    binShapes[proc]=finalNomShape.Clone('%sUp'%systVar)
                    binShapes[proc].SetBinContent(xbin,val+unc)
                    saveToShapesFile(outFile,binShapes,binShapes[proc].GetName())

                    binShapes[proc]=finalNomShape.Clone('%sDown'%systVar)
                    binShapes[proc].SetBinContent(xbin,ROOT.TMath.Max(val-unc,1e-3))
                    saveToShapesFile(outFile,binShapes,binShapes[proc].GetName())

                    #write to datacard
                    datacard.write('%32s shape'%systVar)
                    for sig in mainSignalList:
                        if proc==sig:
                            datacard.write('%15s'%'1')
                        else:
                            datacard.write('%15s'%'-')
                    for sig in altSignalList:
                        if proc==sig:
                            datacard.write('%15s'%'1')
                        else:
                            datacard.write('%15s'%'-')
                    for iproc in exp:
                        if iproc in mainSignalList+altSignalList : continue
                        if iproc==proc:
                            datacard.write('%15s'%'1')
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')

                finalNomShape.Delete()


        #rate systematics: these are fixed values common to all processes
        print '\t rate systematics',len(rateSysts)
        for syst,val,pdf,whiteList,blackList in rateSysts:
            if '*CH*' in syst : syst=syst.replace('*CH*',lfs)
            datacard.write('%32s %8s'%(syst,pdf))
            entryTxt=''
            try:
                entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
            except:
                entryTxt='%15s'%('%3.3f'%val)
            for sig in mainSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                if proc in mainSignalList+altSignalList : continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        #weighting systematics
        print '\t weighting systematics',len(weightingSysts)
        for syst,weightList,whiteList,blackList,shapeTreatment,nsigma in weightingSysts:
            if '*CH*' in syst : syst=syst.replace('*CH*',lfs)

            isGen = any("gen" in twght for twght in weightList)
            isGen = isGen or 'toppt' in weightList

            # jes has annoying formatting different

            #get shapes and adapt them
            iexpUp,iexpDn=None,None
            altExp,altExpUp,altExpDn=None,None,None
            if len(weightList)==1:
                if 'jes' in weightList[0] :
                    jesNum=weightList[0].replace('jes','')
                    _,iexpUp=getDistsForHypoTest(cat,rawSignalList,opt,"","jesup_"+jesNum,isGen)
                    _,iexpDn=getDistsForHypoTest(cat,rawSignalList,opt,"","jesdn_"+jesNum,isGen)
                else :
                    _,iexpUp=getDistsForHypoTest(cat,rawSignalList,opt,"",weightList[0]+"up",isGen)
                    _,iexpDn=getDistsForHypoTest(cat,rawSignalList,opt,"",weightList[0]+"dn",isGen)

                    if syst=='tttoppt':
                        complCat= cat.replace('lowpt','highpt') if 'lowpt' in cat else cat.replace('highpt','lowpt')
                        _,altExp   = getDistsForHypoTest(complCat,rawSignalList,opt)
                        _,altExpUp = getDistsForHypoTest(complCat,rawSignalList,opt,"",weightList[0]+"up",isGen)

                        #reset the down variation to be the nominal one
                        iexpDn=exp
                        altExpDn=altExp

            else:

                #put all the shapes in a 2D histogram
                iexp2D={}
                for iw in xrange(0,len(weightList)):
                    w=weightList[iw]
                    _,kexp=getDistsForHypoTest(cat,rawSignalList,opt,"",w,isGen)
                    for proc in kexp:
                        nbins=kexp[proc].GetNbinsX()
                        if not proc in iexp2D:
                            name =kexp[proc].GetName()+'2D'
                            title=kexp[proc].GetTitle()
                            xmin =kexp[proc].GetXaxis().GetXmin()
                            xmax =kexp[proc].GetXaxis().GetXmax()
                            nReplicas=len(weightList)
                            iexp2D[proc]=ROOT.TH2D(name,title,nbins,xmin,xmax,nReplicas,0,nReplicas)
                            iexp2D[proc].SetDirectory(0)
                        for xbin in xrange(0,nbins+2):
                            iexp2D[proc].SetBinContent(xbin,iw+1,kexp[proc].GetBinContent(xbin))

                #create the up/down variations
                iexpUp,iexpDn={},{}
                for proc in iexp2D:

                    #create the base shape
                    if not proc in iexpUp:
                        tmp=iexp2D[proc].ProjectionX("tmp",1,1)
                        tmp.Reset('ICE')
                        nbinsx=tmp.GetXaxis().GetNbins()
                        xmin=tmp.GetXaxis().GetXmin()
                        xmax=tmp.GetXaxis().GetXmax()
                        iexpUp[proc]=ROOT.TH1F(iexp2D[proc].GetName().replace('2D','up'),proc,nbinsx,xmin,xmax)
                        iexpUp[proc].SetDirectory(0)
                        iexpDn[proc]=ROOT.TH1F(iexp2D[proc].GetName().replace('2D','dn'),proc,nbinsx,xmin,xmax)
                        iexpDn[proc].SetDirectory(0)
                        tmp.Delete()

                    #project each bin shape for the different variations
                    for xbin in xrange(0,iexp2D[proc].GetNbinsX()+2):
                        tmp=iexp2D[proc].ProjectionY("tmp",xbin,xbin)
                        tvals=numpy.zeros(tmp.GetNbinsX())
                        for txbin in xrange(1,tmp.GetNbinsX()+1) : tvals[txbin-1]=tmp.GetBinContent(txbin)

                        #mean and RMS based
                        if 'PDF' in syst:
                            mean=numpy.mean(tvals)
                            rms=numpy.std(tvals)
                            iexpUp[proc].SetBinContent(xbin,mean+rms)
                            iexpDn[proc].SetBinContent(xbin,ROOT.TMath.Max(mean-rms,1.0e-4))

                        #envelope based
                        else:
                            imax=numpy.max(tvals)
                            if iexpUp[proc].GetBinContent(xbin)>0 : imax=ROOT.TMath.Max(iexpUp[proc].GetBinContent(xbin),imax)
                            iexpUp[proc].SetBinContent(xbin,imax)

                            imin=numpy.min(tvals)
                            if iexpDn[proc].GetBinContent(xbin)>0 : imin=ROOT.TMath.Min(iexpDn[proc].GetBinContent(xbin),imin)
                            iexpDn[proc].SetBinContent(xbin,imin)

                        tmp.Delete()


                    #all done, can remove the 2D histo from memory
                    iexp2D[proc].Delete()

            #check the shapes
            iRateVars={}
            if shapeTreatment>0:
                for proc in iexpUp:
                    nbins=iexpUp[proc].GetNbinsX()

                    n=exp[proc].Integral(0,nbins+2)
                    nUp=iexpUp[proc].Integral(0,nbins+2)
                    nDn=iexpDn[proc].Integral(0,nbins+2)

                    #cases where the rate between high/low pt needs to be taken into account
                    upSF,dnSF=1.0,1.0
                    if altExp and altExpUp and altExpDn:
                        ntot   = n+altExp[proc].Integral(0,nbins+2)
                        ntotUp = nUp+altExpUp[proc].Integral(0,nbins+2)
                        ntotDn = nDn+altExpDn[proc].Integral(0,nbins+2)

                        upSF=ntot/ntotUp
                        dnSF=ntot/ntotDn

                    #normalize shapes to nominal expectations
                    if nUp>0: iexpUp[proc].Scale(upSF*n/nUp)
                    if nDn>0: iexpDn[proc].Scale(dnSF*n/nDn)

                    #save a rate systematic from the variation of the yields
                    if n==0 : continue
                    nvarUp=max(nUp/n,nDn/n)
                    nvarDn=min(nUp/n,nDn/n)

                    iRateUnc=buildRateUncertainty(nvarDn,nvarUp)
                    if iRateUnc: iRateVars[proc]=iRateUnc


            #write the shapes to the ROOT file
            saveToShapesFile(outFile,iexpUp,('%s_%s_%sUp'%(cat,opt.dist,syst)),opt.rebin)
            saveToShapesFile(outFile,iexpDn,('%s_%s_%sDown'%(cat,opt.dist,syst)),opt.rebin)

            #fill in the datacard
            datacard.write('%32s %8s'%(syst,'shape'))
            entryTxt='%15s'%('%3.3f'%nsigma)
            for sig in mainSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                if proc in mainSignalList+altSignalList : continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            #write the rate systematics as well
            if shapeTreatment!=2: continue
            if len(iRateVars)==0: continue
            datacard.write('%32s %8s'%(syst+'Rate',pdf))
            for sig in mainSignalList:
                if sig in iRateVars and ((len(whiteList)==0 and not sig in blackList and not '-'+sig in blackList) or sig in whiteList):
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList if opt.useAltRateUncs else mainSignalList:
                if sig in iRateVars and ((len(whiteList)==0 and not sig in blackList and not '-'+sig in blackList) or sig in whiteList):
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                if proc in mainSignalList+altSignalList : continue
                if proc in iRateVars and ((len(whiteList)==0 and not proc in blackList and not '-'+proc in blackList) or proc in whiteList):
                    datacard.write('%15s'%iRateVars[proc])
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')


        #systematics from dedicated samples
        print '\t simulated systematics',len(fileShapeSysts)
        for syst,procsAndSamples,shapeTreatment,nsigma in fileShapeSysts:

            if '*CH*' in syst : syst=syst.replace('*CH*',lfs)

            iexpUp,iexpDn={},{}
            for proc in procsAndSamples:
                samples=procsAndSamples[proc]

                hyposToGet=[(opt.mainHypo if proc in rawSignalList else 100.0)]
                isSignal=False
                if proc in rawSignalList:
                    isSignal=True
                    hyposToGet.append( opt.altHypo )

                jexpDn,jexpUp=None,None
                for hypo in hyposToGet:
                    if len(samples)==2:
                        _,jexpDn=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f'%(cat,opt.dist,hypo),samples[0])
                        _,jexpUp=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f'%(cat,opt.dist,hypo),samples[1])
                    else:
                        _,jexpUp=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f'%(cat,opt.dist,hypo),samples[0])

                    newProc=proc
                    if isSignal:
                        newProc=('%sw%.0f'%(proc,hypo)).replace('.','p')
                    jexpUp.values()[0].SetName(newProc)
                    iexpUp[newProc]=jexpUp.values()[0]

                    #if down variation is not found, mirror it
                    try:
                        jexpDn.values()[0].SetName(newProc)
                        iexpDn[newProc]=jexpDn.values()[0]
                    except:
                        idnHisto=jexpUp.values()[0].Clone()
                        idnHisto.SetDirectory(0)
                        for xbin in xrange(0,idnHisto.GetNbinsX()+2):
                            nomVal=exp[newProc].GetBinContent(xbin)
                            newVal=idnHisto.GetBinContent(xbin)
                            diff=ROOT.TMath.Abs(newVal-nomVal)
                            #if 'tWttInterf' not in syst :
                            #    if newVal>nomVal: nomVal-= ROOT.TMath.Max(diff,1e-4)
                            #    else: nomVal+=diff
                            idnHisto.SetBinContent(xbin,nomVal)
                        iexpDn[newProc]=idnHisto

            #check the shapes
            iRateVars={}
            if shapeTreatment>0:
                for proc in iexpUp:
                    nbins=iexpUp[proc].GetNbinsX()

                    #normalize shapes to nominal expectations
                    n=exp[proc].Integral(0,nbins+2)
                    nUp=iexpUp[proc].Integral(0,nbins+2)
                    if nUp>0: iexpUp[proc].Scale(n/nUp)
                    nDn=iexpDn[proc].Integral(0,nbins+2)
                    if nDn>0: iexpDn[proc].Scale(n/nDn)

                    #save a rate systematic from the variation of the yields
                    if n==0 : continue
                    nvarUp=max(nUp/n,nDn/n) # ROOT.TMath.Abs(1-nUp/n)
                    nvarDn=min(nUp/n,nDn/n) # ROOT.TMath.Abs(1-nDn/n)

                    iRateUnc=buildRateUncertainty(nvarDn,nvarUp)
                    if iRateUnc: iRateVars[proc]=iRateUnc


            #write the shapes to the ROOT file
            saveToShapesFile(outFile,iexpUp,('%s_%s_%sUp'%(cat,opt.dist,syst)),opt.rebin)
            saveToShapesFile(outFile,iexpDn,('%s_%s_%sDown'%(cat,opt.dist,syst)),opt.rebin)

            #fill in the datacard
            datacard.write('%32s %8s'%(syst,'shape'))
            for sig in mainSignalList:
                if sig in iexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[sig.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if sig in iexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[sig.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                if proc in mainSignalList+altSignalList : continue
                if proc in iexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[proc.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            #write the rate systematics as well
            if shapeTreatment!=2: continue
            if len(iRateVars)==0 : continue
            datacard.write('%32s %8s'%(syst+'Rate',pdf))
            for sig in mainSignalList:
                if sig in iRateVars :
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList if opt.useAltRateUncs else mainSignalList:
                if sig in iRateVars :
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                if proc in mainSignalList+altSignalList : continue
                if proc in iRateVars :
                    datacard.write('%15s'%iRateVars[proc])
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        print '\t ended datacard generation'
        datacard.close()

        if opt.doValidation:
            print '\t running validation'
            for proc in rawSignalList:
                newProc=('%sw%.0f'%(proc,opt.mainHypo)).replace('.','p')
                altProc=('%sw%.0f'%(proc,opt.altHypo)).replace('.','p') if proc=='tbart' else ''
                for uncList in [ 'jes,jer,les_*CH*',
                                 'btag,ltag,pu,trig_*CH*',
                                 'ttPSScale,ttMEqcdscale,ttPDF,tttoppt',
                                 'ttGenerator,ttPartonShower',
                                 'mtop',
                                 'tWttInterf,tWQCDScale'
                                 ]:
                    if 'Singletop' in proc and ('ttPS' in uncList or 'ttGen' in uncList) : continue
                    if 'tbart' in proc and 'tWttInterf' in uncList : continue
                    uncList=uncList.replace('*CH*',lfs)
                    plotter=Popen(['python',
                                   '%s/src/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getShapeUncPlots.py'%(os.environ['CMSSW_BASE']),
                                   '-i','%s/shapes.root'%outDir,
                                   '--cats','%s'%cat,
                                   '--obs', '%s'%opt.dist,
                                   '--proc','%s'%newProc,
                                   '-o','%s'%outDir,
                                   '--alt','%s'%altProc,
                                   '--uncs','%s'%uncList],
                                  stdout=PIPE,
                                  stderr=STDOUT)
                    plotter.communicate()

    return outDir,dataCardList

"""
steer the script
"""
def main():

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(      '--combine',            dest='combine',            help='CMSSW_BASE for combine installation',         default=None,        type='string')
    parser.add_option('-i', '--input',              dest='input',              help='input plotter',                               default=None,        type='string')
    parser.add_option(      '--systInput',          dest='systInput',          help='input plotter for systs from alt samples',    default=None,        type='string')
    parser.add_option('-d', '--dist',               dest='dist',               help='distribution',                                default='minmlb',    type='string')
    parser.add_option(      '--nToys',              dest='nToys',              help='toys to through for CLs',                     default=2000,        type=int)
    parser.add_option('--addBinByBin',              dest='addBinByBin', help='add bin-by-bin stat uncertainties @ threshold',      default=-1,            type=float)
    parser.add_option(      '--rebin',              dest='rebin',       help='histogram rebin factor',                             default=0,             type=int)
    parser.add_option(      '--pseudoData',         dest='pseudoData',         help='pseudo data to use (-1=real data)',           default=100,         type=float)
    parser.add_option(      '--useAltRateUncs',     dest='useAltRateUncs',     help='use rate uncertainties specific to alt. hypothesis', default=False,       action='store_true')
    parser.add_option(      '--replaceDYshape',     dest='replaceDYshape',     help='use DY shape from syst file',                 default=False,       action='store_true')
    parser.add_option(      '--doValidation',       dest='doValidation',       help='create validation plots',                     default=False,       action='store_true')
    parser.add_option(      '--rndmPseudoSF',       dest='rndmPseudoSF',       help='multiply pseudodate by random SF?',           default=False,       action='store_true')
    parser.add_option(      '--pseudoDataFromSim',  dest='pseudoDataFromSim',  help='pseudo data from dedicated sample',           default='',          type='string')
    parser.add_option(      '--pseudoDataFromWgt',  dest='pseudoDataFromWgt',  help='pseudo data from weighting',                  default='',          type='string')
    parser.add_option(      '--mainHypo',           dest='mainHypo',  help='main hypothesis',                                      default=100,         type=float)
    parser.add_option(      '--altHypo',            dest='altHypo',   help='alternative hypothesis',                               default=400,         type=float)
    parser.add_option(      '--altHypoFromSim',     dest='altHypoFromSim',   help='alternative hypothesis from dedicated sample',  default='',          type='string')
    parser.add_option('-s', '--signal',             dest='signal',             help='signal (csv)',                                default='tbart,Singletop',  type='string')
    parser.add_option(      '--removeNuisances',    dest='rmvNuisances',       help='nuisance group to remove (csv)',              default='',  type='string')
    parser.add_option(      '--freezeNuisances',    dest='frzNuisances',       help='nuisance group to freeze (csv)',              default='',  type='string')
    parser.add_option('-c', '--cat',                dest='cat',                help='categories (csv)',
                      default='EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt',
                      type='string')
    parser.add_option('-o', '--output',             dest='output',             help='output directory',                            default='datacards', type='string')
    (opt, args) = parser.parse_args()

    outDir,dataCardList=doDataCards(opt,args)
    scriptname=doCombineScript(opt,args,outDir,dataCardList)
    print 'Running statistical analysis'
    runCombine=Popen(['sh',scriptname],stdout=PIPE,stderr=STDOUT)
    runCombine.communicate()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
