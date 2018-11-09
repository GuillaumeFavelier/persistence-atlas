#!/bin/sh

echo "Gaussians:"
time ./gaussians.sh > /dev/null 2>&1

echo "Vortex street:"
time ./vortexStreet.sh > /dev/null 2>&1

echo "Starting vortex:"
time ./startingVortex.sh > /dev/null 2>&1

echo "Isabella:"
time ./isabella.sh > /dev/null 2>&1

echo "Sea Surface Height:"
time ./seaSurfaceHeight.sh > /dev/null 2>&1
