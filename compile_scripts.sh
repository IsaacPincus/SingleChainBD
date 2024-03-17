#!/bin/bash
#Scripts to compile the code, because Makefiles are stupid

## M3Massive parameters as of Feb 13 2020 ##
FC90=mpif90
#FC90NOMPI = gfortran
#FC = gfortran
#CC = gcc

#OPT = -O3 -LFposix
db="-g -fno-omit-frame-pointer -O3 -Wall -fcheck=all"
#db="-g -fno-omit-frame-pointer -O3"
Net1="-lnetcdff"
Net2="-lnetcdf"
#INCLUDES = -I/usr/local/netcdf/4.4.1.1-openmpi-1.10.7-mlx/include
#FFLAGS = $(db) $(OPT) $(INCLUDES) $(Net2) $(Net1)
#FC90FLAGS = $(FFLAGS)
#SYSLIBS = -L/usr/local/netcdf/4.4.1.1-openmpi-1.10.7-mlx/lib -llapack -lblas -lm

## Monarch parameters as of Jun 8 2020 ##
# All the same as M3Massive, only with these syslibs instead:
#SYSLIBS = -L/usr/local/netcdf/4.4.1.1/lib -llapack -lblas -lm

## Local parameters as of Jun 17 2020 ##
FC90="mpifort"
INCLUDES="-L/usr/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial -I/usr/include -lhdf5"
FFLAGS="${db} ${OPT} ${INCLUDES} ${Net2} ${Net1}"
#FC90FLAGS= $(FFLAGS)
#FCFLAGS= $(FFLAGS)
SYSLIBS="-L/usr/lib -llapack -lblas -lm"

FILELIST="modules.f90 DataOutput.f90 utils.f90 gsipc.f90 properties.f90 InitialPos.f90 sensemble.f90"

#actually compile the thing
run="mpifort -o sens ${FILELIST} ${FFLAGS} ${SYSLIBS}"
echo $run
$run


