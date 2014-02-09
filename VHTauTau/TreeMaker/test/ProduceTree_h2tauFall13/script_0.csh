#!/bin/csh
setenv QUI $PWD
setenv VO_CMS_SW_DIR /opt/exp_soft/cms
source $VO_CMS_SW_DIR/cmsset_default.csh
setenv CASA /cmshome/venditti/TauRun2Analysis/CMSSW_6_2_3_patch1/src/VHTauTau/TreeMaker/test/ProduceTree_h2tauFall13
echo $CASA
cd $CASA
setenv SCRAM_ARCH slc5_amd64_gcc472
eval `scramv1 runtime -csh`

cmsRun test_data0.py >& /lustre/cms/store/user/rosma/UpgradePh1_FullSim_16July/GluGluToHToTauTau_Fall13_Tree/0.txt


