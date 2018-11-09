#!/bin/sh

if [ ! -d isabella ];
then
  mkdir isabella
fi
cd isabella
../install/bin/persistenceAtlasCmd -i ../data/isabella_velocity_goodEnsemble.vti -numberOfClusters 3
