#!/bin/bash

echo ===========================
BASE="./AnalysisInput/"
YEAR=2016
PERIOD="LHC16h"
SUBDIR="muon_caloLego"
TRAIN=440
TIME=_20160601-1606

cd $BASE
mkdir $PERIOD
mkdir $PERIOD/$SUBDIR
cd $PERIOD/$SUBDIR
for i in 254381 254394 254395 254396 254476 254479 254604 254606 254607 254608 254621 254629 254630 254632 254640 254644 254646 254648 254649 254651 254652 254653 254654 255009 255010 255011
do

inFILE="alien:///alice/data/$YEAR/$PERIOD/000"${i}"/muon_calo_pass1/PWGPP/PP_EMCAL_Calibration/$TRAIN$TIME/AnalysisResults.root"
echo  $inFILE

outputFILE="file:"${i}".root"
echo $outputFILE

alien_cp $inFILE $outputFILE

done

