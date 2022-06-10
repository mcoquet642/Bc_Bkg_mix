#!/bin/bash
export RAPIDSIM_ROOT=/local/home/mc262512/RapidSim
echo $RAPIDSIM_ROOT

#$RAPIDSIM_ROOT/build/src/RapidSim.exe Dm2K0mumnmu 500000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Dp2K0bmupnmu 500000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Kmmumnu 1000000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Kpmupnu 1000000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Pimmumnu 1000000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Pipmupnu 1000000 1

$RAPIDSIM_ROOT/build/src/RapidSim.exe Bcm2Jpsimunu 16000 1
$RAPIDSIM_ROOT/build/src/RapidSim.exe Bcp2Jpsimunu 16000 1

#$RAPIDSIM_ROOT/build/src/RapidSim.exe Bcm2Jpsimunu_PTCUT 100000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe Bcp2Jpsimunu_PTCUT 100000 1

#$RAPIDSIM_ROOT/build/src/RapidSim.exe B02JpsiKpPim 500000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe B02JpsiKpPim_PTCUT 500000 1

#$RAPIDSIM_ROOT/build/src/RapidSim.exe JpsiMumu_PTCUT 500000 1
#$RAPIDSIM_ROOT/build/src/RapidSim.exe JpsiMumu 500000 1
