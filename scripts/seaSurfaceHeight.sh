#!/bin/sh

if [ ! -d seaSurfaceHeight ];
then
  mkdir seaSurfaceHeight
fi
cd seaSurfaceHeight
../install/bin/persistenceAtlasCmd -i ../data/seaSurfaceHeightGoodEnsemble.vti -spread 0.01 -threshold 0.1 | tee seaSurfaceHeight.log
