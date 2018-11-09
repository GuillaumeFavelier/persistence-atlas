#!/bin/sh

if [ ! -d startingVortex ];
then
  mkdir startingVortex
fi
cd startingVortex
../install/bin/persistenceAtlasCmd -i ../data/startingVortexGoodEnsemble.vti -numberOfNeighbors 4 -critical 2 | tee startingVortex.log
