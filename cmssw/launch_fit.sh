#!/bin/bash

WORKDIR=$PWD

# trick for SWAN: unset previous python env
unset PYTHONPATH
unset PYTHONHOME
# activation of el7 container if needed
#cmssw-el7 -B /eos -B /media/
# load CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
export RELEASE=CMSSW_11_3_4
if [ -r $RELEASE/src ] ; then
  echo found $RELEASE
else
  echo please setup $RELEASE env first
  exit 1
fi
cd $RELEASE/src
eval `scram runtime -sh`

# launch the fit
cd $WORKDIR
python $WORKDIR/cmssw/fit.py $@
