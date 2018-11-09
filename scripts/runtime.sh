#!/bin/sh

pwd

echo "Gaussians:"
time scripts/gaussians.sh > /dev/null 2>&1

echo "Vortex street:"
time scripts/vortexStreet.sh > /dev/null 2>&1

echo "Starting vortex:"
time scripts/startingVortex.sh > /dev/null 2>&1

echo "Isabella:"
time scripts/isabella.sh > /dev/null 2>&1

echo "Sea Surface Height:"
time scripts/seaSurfaceHeight.sh > /dev/null 2>&1
