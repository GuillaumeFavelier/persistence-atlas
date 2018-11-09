#!/bin/sh

if [ ! -d vortexStreet ];
then
  mkdir vortexStreet
fi
cd vortexStreet
../install/bin/persistenceAtlasCmd -i ../data/vortexStreetGoodEnsemble2.vti -spread 0.03 -critical 2 -mcpThreshold 0.015
