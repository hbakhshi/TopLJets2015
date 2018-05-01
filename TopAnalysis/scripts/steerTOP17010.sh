#!/bin/bash

WHAT=$1; 
if [ "$#" -lt 1 ]; then 
    echo "steerTOP17010.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW/HYPOTEST> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        DY           - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo "        HYPOTEST     - create the datacards, steering scripts for hypothesis testing and submit to batch";
    echo "        PLOTHYPOTEST - create summaries of the hypothesis tests";
    echo "        WWWHYPOTEST  - move summaries to afs-web-based area"
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=workday
lumi=35922
lumiSpecs=""
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/b312177
dataeosdir=/store/cmst3/group/top/ReReco2016/be52dbe_03Feb2017
summaryeosdir=/store/cmst3/group/top/TOP-17-010-final-v2/
COMBINERELEASE=${HOME}/CMSSW_8_1_0/src/
outdir=/afs/cern.ch/work/${myletter}/${whoami}/TOP-17-010-final-v2/
anadir=${outdir}/$2
wwwdir=${HOME}/www/TOP-17-010/


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    TEST )
        file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_SingleTbar_tW/MergedMiniEvents_0_ext0.root
        #file=root://eoscms//eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root
	python scripts/runLocalAnalysis.py -i ${file} \
            -q local -o /tmp/`whoami`/MC13TeV_SingleT_tW_test.root --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts;
        python scripts/runTopWidthAnalysis.py -i /tmp/`whoami`/MC13TeV_SingleT_tW_test.root -o /tmp/`whoami`/Chunks -q local;
        ;;
    SEL )
        commonOpts="-q ${queue} -o ${summaryeosdir} --era era2016 -m TOP-17-010::RunTop17010 --ch 0 --runSysts --skipexisting";
	python scripts/runLocalAnalysis.py -i ${eosdir} ${commonOpts}     --only MC --farmappendix TOP17010MC;
	python scripts/runLocalAnalysis.py -i ${dataeosdir} ${commonOpts} --only Data --farmappendix TOP17010Data;
	;;
    MERGESEL )
	mkdir -p ${outdir}
	./scripts/mergeOutputs.py /eos/cms${summaryeosdir} True ${outdir};	
	;;
    PLOTSEL )
        commonOpts="-i ${outdir} --puNormSF puwgtctr  -j data/era2016/samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc} --doDataOverMC"
	python scripts/plotter.py ${commonOpts} 
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    TESTANA )
	python scripts/runTopWidthAnalysis.py -i root://eoscms//eos/cms/${summaryeosdir}/Chunks/MC13TeV_TTJets_12.root -o ${outdir}/analysis/Chunks -q local;
        ;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${outdir}/analysis/Chunks -q ${queue} --only MC13TeV_TTJets;
        #python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${outdir}/analysis/Chunks -q ${queue} --only Data13TeV_Single,Data13TeV_Double --farm TOP17010DataANA;
        #python scripts/runTopWidthAnalysis.py -i ${summaryeosdir}/Chunks -o ${outdir}/analysis/Chunks -q ${queue} --only MuonEG --farm TOP17010DataMuEGANA;
	;;
    CHECKANA )
        for FARM in TOP17010ANA TOP17010DataANA TOP17010DataMuEGANA; do
            python scripts/checkAnalysisIntegrity.py ${CMSSW_BASE}/${FARM} ${outdir}/analysis/Chunks;
        done
        ;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
        opts="-j data/era2016/samples.json  -l ${lumi} ${lumiSpecs} --onlyData --doDataOverMC --mcUnc ${lumiUnc}"
	python scripts/plotter.py -i ${outdir}/analysis ${opts} --only mll -o dy_plotter.root; 
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/ > ${outdir}/analysis/plots/dysf.dat;
        opts="${opts} --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck --doDataOverMC --mcUnc ${lumiUnc}"
	python scripts/plotter.py -i ${outdir}/analysis ${opts} --only mll,evcount,dphilb,drlb,met,njets,ptlb,incmlb_w100 -o dysf_plotter.root; 
	python scripts/plotter.py -i ${outdir}/analysis ${opts} --only count --saveTeX -o count_plotter.root;
	;;
    PLOT )
        opts="-l ${lumi} ${lumiSpecs} --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck --mcUnc ${lumiUnc} --silent"
        python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/samples.json       ${opts};
        #python scripts/plotter.py -i ${outdir}/analysis  -j data/era2016/syst_samples.json  ${opts} -o syst_plotter.root;
        ;;
    COMBPLOT )
	#combined plots
	#python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb       EE1b,EE2b,MM1b,MM2b,EM1b,EM2b     #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb       EE1b,MM1b,EM1b                    #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb       EE2b,MM2b,EM2b                    #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100 EE1blowpt,MM1blowpt,EM1blowpt    #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100 EE2blowpt,MM2blowpt,EM2blowpt    #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100 EE1bhighpt,MM1bhighpt,EM1bhighpt #${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_w100 EE2bhighpt,MM2bhighpt,EM2bhighpt #${outdir}/analysis/plots/plotter.root
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	mkdir -p ${wwwdir}/comb
        cp plots/*.{root,png,pdf} ${wwwdir}/comb 
        cp test/index.php ${wwwdir}/comb
	;;
    HYPOTEST ) 

    mainHypo=(100)

	    #"EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt"
	    #"EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt"
	    #"EM1bhighpt"
        #"EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt"

	CATS=(
        "EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt"
        )
        #"EM1blowpt"
        #"EM1bhighpt"
        #"EM2blowpt"
        #"EM2bhighpt"
        #"MM1blowpt"
        #"MM1bhighpt"
        #"MM2blowpt"
        #"MM2bhighpt"

    TAGS=(
        "inc_scan_ARC_test"
        )
        #"inc_scan_ARC__EM1blowpt"
        #"inc_scan_ARC__EM1bhighpt"
        #"inc_scan_ARC__EM2blowpt"
        #"inc_scan_ARC__EM2bhighpt"
        #"inc_scan_ARC__MM1blowpt"
        #"inc_scan_ARC__MM1bhighpt"
        #"inc_scan_ARC__MM2blowpt"
        #"inc_scan_ARC__MM2bhighpt"

    altHypo=(20 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)        
    #altHypo=(100) 

    altMass=(1710 1712 1714 1716 1718 1720 1722 1724 1725 1726 1728 1730 1732 1734 1736 1738 1740)
    #altMass=(100) 

	#data=(-1 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)
    data=(100)

    expAltHypo=("nan") #"meq166p5" "meq169p5" "meq171p5" "meq173p5" "meq175p5" "meq178p5")

    #nuisanceFreeze=("sel,trig_*CH*,lumi_13TeV,DYnorm_*CH*,Wnorm_th,tWnorm_th,VVnorm_th,tbartVnorm_th,ees,mes,jer,ltag,btag,bfrag,semilep,pu,tttoppt,ttMEqcdscale,ttPDF,jes,st_wid,UE,CR,hdamp,ISR,FSR,mtop,tWttInterf,tWMEScale") 
    #nuisanceFreeze=("sel,trig_*CH*" "lumi_13TeV" "DYnorm_*CH*" "Wnorm_th" 
    #            "tWnorm_th" "VVnorm_th" "tbartVnorm_th" 
    #            "ees" "mes" "jer" "ltag" "btag" "bfrag" "semilep"
    #            "pu" "tttoppt" "ttMEqcdscale" "ttPDF"
    #            "jes" "st_wid" "UE" "CR" 
    #            "hdamp" "ISR" "FSR" "mtop" 
    #            "tWttInterf" "tWMEScale" "all") 
    #nuisanceFreeze=("jes0" "jes1" "jes2" "jes3" "jes4" "jes5" "jes6" "jes7" "jes8" "jes9" "jes10" 
    #                "jes11" "jes12" "jes13" "jes14" "jes15" "jes16" "jes17" "jes18" "jes19" 
    #                "jes20" "jes21" "jes22" "jes23" "jes24" "jes25" "jes26" "jes27" "jes28")
    #nuisanceFreeze=("CR,UE,ISR,FSR,hdamp,pu,ttMEqcdscale,tttoppt")
    nuisanceFreeze=("nan")
    #nuisanceRemove=("nan")

    queue=1nw
    
	for h in ${altHypo[@]}; do
	for mm in ${altMass[@]}; do
    for mh in ${mainHypo[@]}; do
	    for f in ${!nuisanceFreeze[*]}; do
	    for e in ${expAltHypo[@]}; do
	    for d in ${data[@]}; do
		for i in ${!TAGS[*]}; do
                    
            #iNuisRmv=${nuisanceRemove[${f}]}
            iNuisFrz=${nuisanceFreeze[${f}]}
		    icat=${CATS[${i}]}
		    itag=${TAGS[${i}]}
		    
		    cmd="python test/TopWidthAnalysis/runHypoTestDatacards.py"
		    cmd="${cmd} --combine ${COMBINERELEASE}"
		    cmd="${cmd} --mainHypo=${mh} --altHypo ${h} --pseudoData=${d}"
            cmd="${cmd} -s tbart"
            #cmd="${cmd} --rndmPseudoSF"
            if [[ "${e}" != "nan" ]] ; then 
    		    cmd="${cmd} --altHypoFromSim ${e}"		    
            fi
            if [[ "${iNuisFrz}" != "nan" ]] ; then 
    		    cmd="${cmd} --freezeNuisances ${iNuisFrz}"		    
            fi
            #if [[ "${iNuisRmv}" != "nan" ]] ; then 
    		#    cmd="${cmd} --removeNuisances ${iNuisRmv}"		    
            #fi
		    cmd="${cmd} --dist incmlb"		    
		    cmd="${cmd} --tmass 1725"		    
		    cmd="${cmd} --alttmass ${mm}"		    
		    #cmd="${cmd} --dist tmass"		    
		    cmd="${cmd} --nToys 2000"
		    #cmd="${cmd} -i /eos/cms/${summaryeosdir}/plotter/plotter.root"
		    #cmd="${cmd} --systInput /eos/cms/${summaryeosdir}/plotter/syst_plotter.root"
		    cmd="${cmd} -i /afs/cern.ch/work/e/ecoleman/TOP-17-010-final-v2/plotter.root"
		    cmd="${cmd} --systInput /afs/cern.ch/work/e/ecoleman/TOP-17-010-final-v2/syst_plotter.root"
		    #cmd="${cmd} --pseudoDataFromSim=widthx0p5"
		    cmd="${cmd} -c ${icat}"
		    cmd="${cmd} --rebin 2"            
            #cmd="${cmd} --doValidation"
		    #if [ "$h" == "220" ]; then
		    #if [ "$d" == "-1" ]; then
		    #echo "    validation will be included"
		    #cmd="${cmd} --doValidation"
		    #fi
		    #fi
                    
		    echo "Submitting ($mh,$h,$d,$itag,$icat)"		
		    stdcmd="${cmd} -o ${outdir}/datacards_${itag}"

            # add name for frozen syst
            #if [[ "${iNuisFrz}" != "nan" ]] ; then 
    		#    stdcmd="${stdcmd}Frz_${iNuisFrz}"		    
            #fi

            # add name for alternate mainHypo
            #if [[ "${mh}" != "100" ]] ; then 
    		#    stdcmd="${stdcmd}_${mh}"		    
            #fi

            # submit job
    		stdcmd="${stdcmd}/"		    
            echo ${stdcmd}
    	    bsub -q ${queue} sh ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${stdcmd};
            sleep 1 
           
		    #if [ "$itag" == "inc" ]; then
			#if [ "$d" == "100" ]; then
			#    #echo "    injecting pseudo-data from nloproddec"
			#    #nlocmd="${cmd} --pseudoDataFromWgt nloproddec -o ${outdir}/datacards_${itag}_nloproddec"
			#    #bsub -q ${queue} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${nlocmd};   
			#    #echo "    injecting pseudo-data from widthx4"
			#    #width4cmd="${cmd} --pseudoDataFromSim=t#bar{t}_widthx4 -o ${outdir}/datacards_${itag}_widthx4"
			#    #bsub -q ${queue} sh ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh ${width4cmd};
			#fi
		    #fi
		done
        done
        done
        done
        done
        done
	    done
	;;
    PLOTHYPOTEST )
	summaryScript=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/summarizeHypoTestResults.py;
	cd ${COMBINERELEASE}/
	eval `scramv1 r -sh`
	cd -
	TAGS=("inc") # "ll" "em" "inc_nloproddec" "inc_widthx4")
	for i in ${!TAGS[*]}; do
	    itag=${TAGS[${i}]}
	    python ${summaryScript} -i ${outdir}/datacards_${itag}/ --doCLs --recreateCLsSummary --doNuisances --doFitSummary	
	done	
	;;
    WWWHYPOTEST )
	TAGS=("_ll" "_em" "_inc" "_inc_nloproddec" "_inc_widthx4")
        for i in ${!TAGS[@]}; do
	    itag=${TAGS[${i}]}
            mkdir -p ${wwwdir}/hypo${itag}
            cp ${outdir}/datacards${itag}/*.{png,pdf} ${wwwdir}/hypo${itag};
            cp ${outdir}/datacards${itag}/hypotest_100vs220_data/*.{png,pdf} ${wwwdir}/hypo${itag};
            cp test/index.php ${wwwdir}/hypo${itag}
	done
	;;
    IMPACTS )
        echo "Starting nuisance impact estimation..."
        cd /afs/cern.ch/work/e/ecoleman/CMSSW_8_1_0/
        eval `scramv1 r -sh`
        cd ${outdir}/$2/

        combineTool.py -M Impacts -d workspace.root -m 172.5 --doInitialFit --robustFit 1
        combineTool.py -M Impacts -d workspace.root -m 172.5 --robustFit 1 --doFits --job-mode lxbatch --task-name nuis-impacts --sub-opts='-q 1nw'
        ;;
    FINALIMP )
        echo "Finalizing nuisance impacts..."
        cd /afs/cern.ch/work/e/ecoleman/CMSSW_8_1_0/
        eval `scramv1 r -sh`
        cd ${outdir}/$2/
       
        combineTool.py -M Impacts -d workspace.root -m 172.5 -o impacts.json
        plotImpacts.py -i impacts.json -o impacts
        ;;
esac
