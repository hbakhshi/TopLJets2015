import ROOT
import os
import sys
import argparse
import pickle
import copy

#sigma=1pb distributed accross crossing angles 
#NB this does not sum to 1 as we don't use all crossing angles in the analysis
SIGNALXSECS       = {120:0.269,130:0.273,140:0.143,150:0.293}
PHOTONSIGNALXSECS = {120:0.372,130:0.295,140:0.162,150:0.171}
VALIDLHCXANGLES   = [120,130,140,150]
CH_DICT           = {'169':'zmm','121':'zee','22':'g'}
CH_TITLE_DICT     = {'169':'Z#rightarrow#mu#mu','121':'Z#rightarrowee','22':'#gamma'}


def smoothMissingMass(h,smoothRanges=[#(-300,   300,  'pol2'),
                                      ( 1500,  2500, 'gaus'),
                                      (-2500, -1500, 'gaus')]):

    """smooths missing mass spectrum in specified ranges"""

    smoothH=h.Clone(h.GetName()+'_smooth')

    hxmin=h.GetXaxis().GetXmin()
    hxmax=h.GetXaxis().GetXmax()

    for xran in smoothRanges:
        
        xmin,xmax,fname=xran
        xmin=min(hxmax,max(hxmin,xmin))
        xmax=min(hxmax,max(hxmin,xmax))

        if xmin==xmax :
            print 'Skipping', xran,'for',h.GetName()
            continue

        imin,imax=h.GetXaxis().FindBin(xmin),h.GetXaxis().FindBin(xmax)
        h.Fit(fname,'RQM+','',xmin,xmax)
        ffunc=h.GetFunction(fname)
        try:        
            for xbin in range(imin,imax) :
                xcen=h.GetXaxis().GetBinCenter(xbin)
                smoothH.SetBinContent(xbin,ffunc.Eval(xcen))
            ffunc.Delete()
        except:
            print '[Warn] no smoothing applied in',xmin,xmax,'to',h.GetName()

    return smoothH


def defineProcessTemplates(histos,norm=False):

    """defines the nominal template and the variations and checks fo 0's in the histograms"""

    templates=[]

    #nominal
    templates.append( histos[0] )
    nomStats=templates[-1].Integral()

    #systematic variations
    #if Up/Down already in the name store directly updating the name
    #if not, mirror the variation given 
    for i in xrange(1,len(histos)):        

        templates.append( histos[i] )
        if norm : templates[-1].Scale(nomStats/histos[i].Integral())
        
        key=templates[-1].GetName()
        if not 'Up' in key and not 'Down' in key :

            templates[-1].SetName(key+'Up')
            templates.append( histos[0].Clone(key+'Down') )
            
            #diff
            diffH=templates[-1].Clone('ratio')
            diffH.Add(templates[-2],-1)
            for xbin in range(templates[0].GetNbinsX()):
                diffVal=diffH.GetBinContent(xbin+1)
                templates[-1].SetBinContent(xbin+1,max(1e-6,diffVal+histos[0].GetBinContent(xbin+1)))                
            diffH.Delete()

            if norm : templates[-1].Scale(nomStats/histos[i].Integral())
    
    #don't leave bins with 0's
    for h in templates:
        h.SetDirectory(0)
        iStats=h.Integral()
        if iStats>0: h.Scale(nomStats/iStats)
        for xbin in range(h.GetNbinsX()):
            if h.GetBinContent(xbin+1)>0: continue
            h.SetBinContent(xbin+1,1e-6)
            h.SetBinError(xbin+1,1e-6)
            
    return templates

        
def fillBackgroundTemplates(opt):

    """fills the background and observed data histograms"""

    totalBkg={}
    templates=[]
    data_templates=[]

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(opt.input,x) for x in os.listdir(opt.input) if 'Data13TeV' in x]:
        if 'MuonEG' in f : continue
        if opt.chTag.find('zmm')==0:
            if 'Photon' in f or 'DoubleEG' in f : 
                continue
        if opt.chTag.find('zee')==0:
            if 'Muon' in f or 'Photon' in f : 
                continue
        if opt.chTag.find('g_')==0:
            if not 'Photon' in f : 
                continue
        data.AddFile(f)

    #define final preselection cuts
    cuts=opt.presel

    #loop over categories build templates
    for icat in range(len(opt.categs)):

        #apply category cuts
        categCut=opt.categs[icat]
        categCut=cuts if len(categCut)==0 else '%s && %s'%(categCut,cuts)
        catName='%s_%d'%(opt.chTag,icat)

        print '\t\t',catName,categCut

        #background modelling histos
        histos=[]
        data_obs=None
        histodef='mmiss >> h(%d,%f,%f)'%(opt.nbins,opt.mMin,opt.mMax)
        if opt.signed:
            #histodef='(ypp>=0 ? mmiss : -mmiss) >> h(%d,%f,%f)'%(2*opt.nbins,-opt.mMax,opt.mMax)
            histodef='(bosoneta>=0 ? mmiss : -mmiss) >> h(%d,%f,%f)'%(2*opt.nbins,-opt.mMax,opt.mMax)
        data.Draw(histodef,
                  '{0} && mmiss>0 && mixType==0'.format(categCut),
                  'goff')
        h=data.GetHistogram()
        totalBkg[icat]=h.Integral()
        if opt.unblind:
            data_obs=h.Clone('data_obs_'+catName)
            data_obs.SetDirectory(0)
        h.Reset('ICE')

        for name,mixType,pfix in [('bkg_'+catName,                          1, ''),
                                  ('bkg_%s_bkgShape'%catName,               1, 'syst'),
                                  ('bkg_%s_bkgShapeSingleDiffUp'%catName,   2, ''),
                                  ('bkg_%s_bkgShapeSingleDiffDown'%catName, 2, 'syst'),
                              ]:

            templCuts=categCut.replace('csi1',pfix+'csi1')
            templCuts=templCuts.replace('csi2',pfix+'csi2')
            if pfix=='syst':
                templCuts=templCuts.replace('protonCat','systprotonCat')

            histodef='%smmiss >> h(%d,%f,%f)'%(pfix,opt.nbins,opt.mMin,opt.mMax)           
            if opt.signed:
                #histodef='(%sypp>=0 ? %smmiss : -%smmiss) >> h(%d,%f,%f)'%(pfix,pfix,pfix,2*opt.nbins,-opt.mMax,opt.mMax)
                histodef='(bosoneta>=0 ? %smmiss : -%smmiss) >> h(%d,%f,%f)'%(pfix,pfix,2*opt.nbins,-opt.mMax,opt.mMax)
            
            data.Draw(histodef,
                      'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(templCuts,pfix,mixType),
                      'goff')
            h=data.GetHistogram()
            h.Scale(totalBkg[icat]/h.Integral())
            histos.append(h.Clone(name))
            histos[-1].SetDirectory(0)

            #use first histogram in a category as pseudo-data in case we're blinded
            if len(histos)==1 and not opt.unblind :
                data_obs=h.Clone('data_obs_'+catName)
                data_obs.SetDirectory(0)

            h.Reset('ICE')

        #finalize templates
        templates += defineProcessTemplates(histos)
        data_templates.append(data_obs)

    print '\t total background:',totalBkg
    return totalBkg,templates,data_templates


def fillSignalTemplates(mass,signalFile,xsec,opt,fiducialCuts='gencsi1>0.02 && gencsi1<0.16 && gencsi2>0.03 && gencsi2<0.18',postfix=''):

    """fills the signal histograms"""

    totalSig={'fid':{},'outfid':{}}
    templates={'fid':[],'outfid':[]}
    nom_templates={'fid':[],'outfid':[]}

    #import signal events
    data=ROOT.TChain('data')
    data.AddFile(os.path.join(opt.input,signalFile))
    dataWgt=14586.4464/41529.3 #preTS2/total

    dataAlt=ROOT.TChain('data')
    dataAlt.AddFile(os.path.join(opt.input,signalFile).replace('preTS2','postTS2'))

    #common weight exppression
    wgtExpr='ppsEff*wgt*{xsec}*{lumi}'.format(xsec=xsec,lumi=opt.lumi)

    #loop over categories build templates
    for icat in range(len(opt.categs)):

        #apply category cuts
        catName='%s_%d%s'%(opt.chTag,icat,postfix)
        categCut=opt.presel
        if len(opt.categs[icat]): categCut += ' && ' + opt.categs[icat]
        print '\t\t',catName,categCut
               
        #define final preselection cuts and repeat for fiducial/non-fiducial regions
        for sigType in totalSig.keys():
                        
            #signal modelling histograms
            histos=[]
            for name,mixType,pfix,addWgt in [('sig_%s_m%s'%(catName,mass),               1, '',     None),
                                             ('sig_%s_m%s_sigShape'%(catName,mass),      1, 'syst', None),                                             
                                             ('sig_%s_m%s_sigCalibUp'%(catName,mass),    1, '',     None),
                                             ('sig_%s_m%s_sigCalibDown'%(catName,mass),  1, '',     None),
                                             ('sig_%s_m%s_sigPPSEff'%(catName,mass),     1, '',     '(ppsEff>0 ? 1.0+0.04/ppsEff : 1.0)'),
                                             ('sig_%s_m%s_sigPzModel'%(catName,mass),    1, '',     'gen_pzwgtUp')]:

                name=sigType+name
                templCuts=categCut.replace('csi1',pfix+'csi1')
                templCuts=templCuts.replace('csi2',pfix+'csi2')
                if sigType=='outfid':
                    templCuts+=' && !(%s)'%fiducialCuts
                else:
                    templCuts+=' && %s'%fiducialCuts
                if pfix=='syst':
                    templCuts=templCuts.replace('protonCat','systprotonCat')

                histodef='%smmiss >> h(%d,%f,%f)'%(pfix,opt.nbins,opt.mMin,opt.mMax)           
                if opt.signed:
                    #histodef='(%sypp>=0 ? %smmiss : -%smmiss) >> h(%d,%f,%f)'%(pfix,pfix,pfix,2*opt.nbins,-opt.mMax,opt.mMax)
                    histodef='(bosoneta>=0 ? %smmiss : -%smmiss) >> h(%d,%f,%f)'%(pfix,pfix,2*opt.nbins,-opt.mMax,opt.mMax)

                shiftDataWgt=1.0
                if 'sigCalibUp' in name:   shiftDataWgt *=1.03
                if 'sigCalibDown' in name: shiftDataWgt *=0.97
                

                #sum up contributions
                data.Draw(histodef,
                          '{0}*{1}*{2}*({3} && mixType=={4} && {5}mmiss>0)'.format(wgtExpr,
                                                                                   addWgt if addWgt else '1',
                                                                                   dataWgt*shiftDataWgt,
                                                                                   templCuts,
                                                                                   mixType,
                                                                                   pfix),
                          'goff')
                h=data.GetHistogram()

                dataAlt.Draw(histodef.replace('h(','halt('),
                             '{0}*{1}*{2}*({3} && mixType=={4} && {5}mmiss>0)'.format(wgtExpr,
                                                                                      addWgt if addWgt else '1',
                                                                                      (1-dataWgt*shiftDataWgt),
                                                                                      templCuts,
                                                                                      mixType,
                                                                                      pfix),
                             'goff')
                h.Add(dataAlt.GetHistogram())

                histos.append( h.Clone(name) )         
                histos[-1].SetDirectory(0)

                if len(histos)==1:
                    totalSig[sigType][icat]=h.Integral()
                    nom_templates[sigType].append(h.Clone(name+'_sigforpseudodata'))
                    nom_templates[sigType][-1].SetDirectory(0)

                h.Reset('ICE')
            templates[sigType] += defineProcessTemplates(histos)
    
    print '\t total signal:',totalSig
    return totalSig,templates,nom_templates


def writeDataCards(opt,shapesURL):

    """writes the datacard and the workspace"""

    finalState='mm'
    if opt.chTag.find('zee')==0 : finalState='ee'
    if opt.chTag.find('g')==0:    finalState='g'

    #create a card per category
    dcList=[]
    for icat in range(len(opt.categs)):
        cat='%s_%d'%(opt.chTag,icat)
        dcTxt='%s/shapes-parametric.datacard_%s.dat'%(opt.output,cat)
        dcList.append(dcTxt)
        with open(dcTxt,'w') as dc:
            dc.write('#\n')
            dc.write('# datacard was automatically generated with generateWorkspace.py\n')
            dc.write('# the options passed are printed below\n')
            dc.write('# %s\n'%opt)
            dc.write('#\n')
            dc.write('imax *\n')
            dc.write('jmax *\n')
            dc.write('kmax *\n')
            dc.write('-'*50+'\n')
            dc.write('shapes fidsig    * {0} $PROCESS_{1}_m$MASS $PROCESS_{1}_m$MASS_$SYSTEMATIC\n'.format(shapesURL,cat))
            dc.write('shapes outfidsig * {0} $PROCESS_{1}_m$MASS $PROCESS_{1}_m$MASS_$SYSTEMATIC\n'.format(shapesURL,cat))
            dc.write('shapes bkg       * {0} $PROCESS_{1}        $PROCESS_$SYSTEMATIC\n'.format(shapesURL,cat)) 
            dc.write('shapes data_obs  * {0} $PROCESS_{1}\n'.format(shapesURL,cat))
            dc.write('-'*50+'\n')
            dc.write('bin %s\n'%cat)
            dc.write('observation -1\n')
            dc.write('-'*50+'\n')
            dc.write('%15s %15s %15s %15s\n'%('bin',cat,cat,cat))
            dc.write('%15s %15s %15s %15s\n'%('process','fidsig','outfidsig','bkg'))
            dc.write('%15s %15s %15s %15s\n'%('process','0',     '1',        '2'))
            dc.write('%15s %15s %15s %15s\n'%('rate',   '-1',    '-1',       '-1'))
            dc.write('-'*50+'\n')
            
            #uncertainties
            dc.write('lumi                   %8s %15s %15s %15s\n'%('lnN',                '1.027', '-',  '-'))
            dc.write('eff_%s                 %8s %15s %15s %15s\n'%(finalState,  'lnN',   '1.03',  '-',  '-'))
            dc.write('sigShape               %8s %15s %15s %15s\n'%('shape',              '1',     '1',  '-'))           
            dc.write('sigCalib               %8s %15s %15s %15s\n'%('shape',              '1',     '1',  '-'))
            dc.write('sigPPSEff              %8s %15s %15s %15s\n'%('shape',              '1',     '1',  '-'))
            dc.write('sigPzModel             %8s %15s %15s %15s\n'%('shape',              '1',     '1',  '-'))
            dc.write('%s_bkgShape            %8s %15s %15s %15s\n'%(cat,'shape',          '-',     '-',  '1')) #uncorrelate background shapes
            dc.write('%s_bkgShapeSingleDiff  %8s %15s %15s %15s\n'%(cat,'shape',          '-',     '-',  '1'))

            #template statistics 
            # threshold=effecting entries (neff~1600 is approximately 2% relative unc.
            # but it's safer to leave it to 0 as this is a bump hunt and any stat effect matters...)
            # include signal=0 (1 would count also the signal entries but it has unknown normalization at start in a search)
            dc.write('{0} autoMCStats 0.0 0\n'.format(cat))
        
            #float the background normalization as well as the signal
            dc.write('mu_bkg_%s       rateParam * bkg       1\n'%cat)
            dc.write('mu_outfidsig rateParam * outfidsig 0\n')

    print '\tshapes available @',shapesURL
    print '\tgenerated the following datacards',dcList

def datacardTask(args):

    #finalize configuration of this specific task
    ch,xangle,opt=args


    setattr(opt,'chTag','%s_a%d'%(CH_DICT[ch],xangle))
    setattr(opt,'chTitle',CH_TITLE_DICT[ch])
    if xangle in VALIDLHCXANGLES:
        setattr(opt,'presel','cat==%s && xangle==%d && %s'%(ch,xangle,opt.preselZ))
    else:
        setattr(opt,'presel','cat==%s && %s'%(ch,opt.preselZ))

    boson='Z'
    if ch=='22':
        boson='gamma'
        if xangle in VALIDLHCXANGLES:
            setattr(opt,'presel','cat==%s && xangle==%d && %s'%(ch,xangle,opt.preselGamma))        
        else:
            setattr(opt,'presel','cat==%s && %s'%(ch,opt.preselGamma))

    #start the output
    outName='shapes_%s_a%d.root'%(ch,xangle)
    shapesURL=os.path.join(opt.output,outName)
    fOut=ROOT.TFile.Open(shapesURL,'RECREATE')

    #define background templates
    print '\t filling background templates and observed data'
    bkgExp,bkgTemplates,data_templates=fillBackgroundTemplates(opt)    
    for h in bkgTemplates:        
        h.SetDirectory(fOut)
        h.Write()

    #parametrize the signal
    print '\t filling signal templates' 
    for m in opt.massList:
        try:
            if xangle in VALIDLHCXANGLES:
                signalFile=opt.sig.format(boson=boson,xangle=xangle,mass=m)
                sigExp,sigTemplates,sigNomTemplates=fillSignalTemplates(mass=m,
                                                                        signalFile=signalFile,
                                                                        xsec=SIGNALXSECS[xangle] if boson=='Z' else PHOTONSIGNALXSECS[xangle],
                                                                        opt=opt)
            else:
                print '\t summing up all crossing angles'
                sigExp,sigTemplates,sigNomTemplates=None,None,None
                for ixangle in VALIDLHCXANGLES:
                    print '\t @%durad'%ixangle
                    signalFile=opt.sig.format(boson=boson,xangle=ixangle,mass=m)
                    if sigExp is None:
                        sigExp,sigTemplates,sigNomTemplates=fillSignalTemplates(mass=m,
                                                                                signalFile=signalFile,
                                                                                xsec=SIGNALXSECS[ixangle] if boson=='Z' else PHOTONSIGNALXSECS[ixangle],
                                                                                opt=opt)
                    else:
                        i_sigExp,i_sigTemplates,i_sigNomTemplates=fillSignalTemplates(mass=m,
                                                                                      signalFile=signalFile,
                                                                                      xsec=SIGNALXSECS[ixangle] if boson=='Z' else PHOTONSIGNALXSECS[ixangle],
                                                                                      opt=opt,
                                                                                      postfix='_add')
                        for k in sigTemplates:
                            for kk in sigExp[k]:
                                sigExp[k][kk] += i_sigExp[k][kk]
                            for i in range(len(sigTemplates[k])):
                                sigTemplates[k][i].Add(i_sigTemplates[k][i])
                                i_sigTemplates[k][i].Delete()
                            for i in range(len(sigNomTemplates[k])):
                                sigNomTemplates[k][i].Add(i_sigNomTemplates[k][i])
                                i_sigNomTemplates[k][i].Delete()

                print '\t Final total signal:',sigExp

            for key in sigTemplates:

                for h in sigTemplates[key]:
                    h.SetDirectory(fOut)
                    h.Write()

                #if blinded add signal to pseudo-data 
                if not opt.unblind and opt.injectMass==m:
                    for i in range(len(sigNomTemplates[key])):
                        data_templates[i].Add(sigNomTemplates[key][i],opt.injectMu)
                        sigNomTemplates[key][i].Delete() #no longer needed

        except Exception as e:
            print '-'*50
            print 'Failed to generate templates at mass=',m,'for boson=',boson
            print e
            print '-'*50

    #now write the data
    for h in data_templates:
        h.SetDirectory(fOut)
        h.Write()

    #all done
    fOut.Close()

    #write summary in datacards
    print '\t writing datacard'
    writeDataCards(opt,os.path.basename(shapesURL))
    

def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05/',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--sig',
                        dest='sig',
                        default='{boson}_m_X_{mass}_xangle_{xangle}_2017_preTS2.root',
                        help='signal point [%default]')
    parser.add_argument('--massList',
                        dest='massList',
                        default='600,660,720,780,800,840,900,960,1000,1020,1080,1140,1200,1260,1320,1380,1400,1440,1500,1560,1600',
                        help='signal mass list (CSV) [%default]')
    parser.add_argument('--injectMass',
                        dest='injectMass',
                        default=None,
                        help='mass to inject in pseudo-data [%default]')
    parser.add_argument('--injectMu',
                        dest='injectMu',
                        default=1.0,
                        help='signal strength of the mass to inject in pseudo-data [%default]')
    parser.add_argument('--preselZ',
                        dest='preselZ', 
                        default='l1pt>30 && l2pt>20 && bosonpt>50',
                        help='preselection for Z categories [default: %default]')
    parser.add_argument('--preselGamma',
                        dest='preselGamma', 
                        default='bosonpt>95',
                        help='preselection for photon categories [default: %default]')
    parser.add_argument('--categs',
                        dest='categs',
                        default='nvtx<20,nvtx>=20',
                        help='Sub-categories [default: %default]')
    parser.add_argument('--signed',
                        dest='signed',
                        default=False,
                        help='use signed missing mass [default %default]',
                        action='store_true')
    parser.add_argument('--lumi',
                        dest='lumi',
                        default=37500.,
                        type=float,
                        help='integrated luminosity [default: %default]')
    parser.add_argument('--mBin',
                        dest='mBin',
                        default=40.,
                        type=float,
                        help='mass bin width [default: %default]')
    parser.add_argument('--mMin',
                        dest='mMin',
                        default=0,
                        type=float,
                        help='minimum missing mass [default: %default]')
    parser.add_argument('--mMax',
                        dest='mMax',
                        default=2000,
                        type=float,
                        help='maximum missing mass [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis/stat',
                        help='Output directory [default: %default]')
    parser.add_argument('--xangles',
                        dest='xangles', 
                        default='120,130,140,150',
                        help='Scan these cross angles (0=inclusive) [default: %default]')
    parser.add_argument('--unblind',
                        dest='unblind', 
                        default=False,
                        action='store_true',
                        help='Use non-mixed data in the final fit [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #configuration
    opt.categs=opt.categs.split(',')
    if len(opt.categs)==0 : opt.categs=[]
    opt.massList=opt.massList.split(',')
    setattr(opt,'nbins', ROOT.TMath.FloorNint( (opt.mMax-opt.mMin)/(opt.mBin)) )

    print '[generateWorkspace]'
    print '\t signal for masses=',opt.massList
    print '\t will apply the following preselection:'
    print '\t\t Z:',opt.preselZ
    print '\t\t gamma:',opt.preselGamma
    print '\t histograms defined as (%d,%f,%f)'%(opt.nbins,opt.mMin,opt.mMax)
    if len(opt.categs) : 
        print '\t will categorize in:',opt.categs
    if opt.unblind:
        print '\t Analysis will be unblinded...'
        print '\t ...how sure are you of what you did?'
    else:
        print '\t Pseudo-data built from background expectations'
        print '\t Signal injected from mass=',opt.injectMass,'with mu=',opt.injectMu

    #prepare output
    os.system('mkdir -p %s'%opt.output)

    task_list=[]
    for ch in CH_DICT.keys():
        for angle in [int(x) for x in opt.xangles.split(',')]: 
            task_list.append( (ch,angle,copy.deepcopy(opt)) )

    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map( datacardTask, task_list )

    print '\t all done, output can be found in',opt.output

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



