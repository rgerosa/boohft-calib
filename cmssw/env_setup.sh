#!/bin/bash

WORKDIR=$PWD

# trick for SWAN: unset previous python env
unset PYTHONPATH
unset PYTHONHOME

#cmssw-el7 -B /eos -B /media/
source /cvmfs/cms.cern.ch/cmsset_default.sh
export RELEASE=CMSSW_11_3_4

if [ -r $RELEASE/src ] ; then
    echo release $RELEASE already exists
else
    scram p CMSSW $RELEASE
    cd $RELEASE/src
    eval `scram runtime -sh`

    ## Install HiggsAnalysis
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cp $WORKDIR/cmssw/data/TagAndProbeExtendedV2.py HiggsAnalysis/CombinedLimit/python/  # copy the model we will use in fit
    cd HiggsAnalysis/CombinedLimit
    git checkout v9.2.1 # recommended tag
    cd ../..
    
    ## Install CombineHarvester
    git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
    cp $WORKDIR/cmssw/data/plot1DScanWithOutput.py CombineHarvester/CombineTools/scripts/
    scram b -j8

    cd ../..
fi
