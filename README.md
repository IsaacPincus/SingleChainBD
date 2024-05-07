Fortran program to simulate a single Brownian polymer chain in shear flow with excluded volume and hydrodynamic interactions.

```gen_input.py``` has an overview of possible code options. To run it, one will need to update the file locations manually.

To run the code on wsl or linux, the following packages are necessary:

```
sudo apt-get install libnetcdf-dev libnetcdff-dev
sudo apt-get install libopenmpi-dev
sudo apt-get install build-essential
sudo apt-get install libhdf5-serial-dev
```

Depending on how your system is set up, you may also need to install lapack and blas libraries, for example:
```
sudo apt-get install libblas-dev liblapack-dev
```

You can compile and run the test cases with:
```
Source compile_tests.sh
./tests
```

There may be one failure, which is related to an experimental feature and can be ignored. You may also have to run compile_tests.sh twice, to first generate the required module files.

Using:
```
Source compile_scripts.sh
```

Should give a single file ```sens``` which can be run using ```mpirun```.

The Matlab script ```Calculate_averages_from_trajectories.m``` can be used to post-process the generated netcdf files and extract bead positions, forces, stresses etc.

Please contact ```ipincus@mit.edu``` or ```isaac.pincus@gmail.com``` if you have questions about the code or if it is not compiling/running correctly.


