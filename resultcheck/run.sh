module load x86/netcdf/4.6.1-with-parallel
module load x86/Intel_Development_Tool/2017
bsub -I -q q_x86_expr -n 4 ./test correctresult/test01_G_T62_t12.pop.r.0001-02-01-00000.nc ../test.nc 
