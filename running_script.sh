#$ bin/bash

rm output.txt
gfortran wavefunction_maha.f90 -o output.o 
#/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/nvfortran -acc wavefunction.f90 -Minfo=accel -gpu=cc80
./output.o >> output.txt
python3 data_plot.py 