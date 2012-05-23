#!/usr/bin/env bash

cd src

gfortran -g -fbounds-check -O3 -c rdn_utils.f90
gfortran -g -fbounds-check -O3 -c core_utils.f90 lattice_utils.f90
gfortran -g -fbounds-check -O3 -o HumbleSimulation core.f90 core_utils.o lattice_utils.o rdn_utils.o

mv -u HumbleSimulation ../
mv -u *.mod ../lib
mv -u *.o ../bin

reset
