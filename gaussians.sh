#!/bin/sh

if [ ! -d gaussians ];
then
  mkdir gaussians
fi
cd gaussians
../install/bin/persistenceAtlasCmd -i ../data/main_example.vtu -numberOfNeighbors 4
