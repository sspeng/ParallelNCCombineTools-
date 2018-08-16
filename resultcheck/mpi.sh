module load x86/netcdf/4.6.1-with-parallel
module load x86/Intel_Development_Tool/2017
mpiicc -fopenmp  -o test $1 -lnetcdf
