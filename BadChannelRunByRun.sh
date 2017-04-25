#!/bin/bash
echo "Run by run analysis"
PERIOD="LHC15o"
TRAIN=758
CODEDIR="Run_BadChannel.C"
BASEDIR=$PWD
BASEFOLDER="AnalysisInput"
BASE=$BASEDIR/$BASEFOLDER
SUFFILE="_AnyINTnoBCFiltered.root"
#SUFFILE="_INT7NoBCFiltered.root"

for j in 244918 244975 244980 244982 244983 245064 245066 245068 245145 245146 245151 245152 245231 245232 245259 245343 245345 245346 245347 245349 245353 245396 245397 245401 245407 245409 245411 245439 245441 245446 245454 245496 245497 245501 245504 245505 245507 245535 245540 245542 245543 245544 245545 245554 245683 245700 245702 245705 245738 245829 245831 245833 245949 245952 245954 245963 246001 246003 246037 246042 246052 246053 246087 246089 246113 246115 246217 246222 246225 246271 246272 246390 246391 246392 246424 246434 246487 246488 246493 246495 246750 246751 246757 246758 246759 246760 246765 246766 246804 246805 246807 246808 246809 246810 246844 246845 246846 246928 246945
do
#for j in 244918 244975 244980 244982 244983 245064 245066 245068 245145 245146 245151 245152 245231 245232 245259
#do

   root -b -q $CODEDIR+\(0,\"$PERIOD\",\"Train_$TRAIN\",\"AnyINTnoBC\",$j,\"$j$SUFFILE\",\"\",\".\"\)
done
