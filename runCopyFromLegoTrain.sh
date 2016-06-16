#!/bin/bash

echo ===========================
BASE="./"
#"/data/Work/EMCAL/BadChannels/"
YEAR=2016
PERIOD="LHC16h"
SUBDIR="muon_caloLego"
TRAIN=440
TIME=_20160601-1606
#LHC16g: "_20160601-1606"
#254124 254126 254128 254147 254148 254149 254174 254196 254199 254204 254205 254293 254302 254303 254304 254330 254331 254332

cd $BASE
mkdir $PERIOD
mkdir $PERIOD/$SUBDIR
cd $PERIOD/$SUBDIR
for i in 254381 254394 254395 254396 254476 254479 254604 254606 254607 254608 254621 254629 254630 254632 254640 254644 254646 254648 254649 254651 254652 254653 254654 255009 255010 255011




  

	
do

#inFILE="alien:/alice/data/2010/LHC10h/000"${i}"/ESDs/pass2/QA52/QAresults.root"
#inFILE="alien:/alice/sim/2013/LHC13d3_plus/"${i}"/QAresults.root"
inFILE="alien:///alice/data/$YEAR/$PERIOD/000"${i}"/muon_calo_pass1/PWGPP/PP_EMCAL_Calibration/$TRAIN$TIME/AnalysisResults.root"

echo  $inFILE
#outputFILE=${i}"passcalo.xml"

#outputFILE="file:/scratch/alicehp2/germain/analysis/AlexisAnalysis/data2011/LegoTrainOutput/LHC11dESD/"${i}".root"
outputFILE="file:"${i}".root"

echo $outputFILE
alien_cp $inFILE $outputFILE

done

