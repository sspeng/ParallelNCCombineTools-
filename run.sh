module load x86/netcdf/4.6.1-with-parallel
module load x86/Intel_Development_Tool/2017
#bsub -I -q q_x86_expr -n 27 ./test ../pop/test01_G_T62_t12.pop.r.0001-02-01-00000.nc
#bsub -I -fs_proj online1 -q q_x86_expr -n 14 ./test ../../testdata/pop/test01_G_T62_t12.pop.r.0001-02-01-00000.nc test.nc
#bsub -I -fs_proj online1 -q q_x86_expr -N 2 -np 7 ./test ../testdata/pop/test01_G_T62_t12.pop.r.0001-02-01-00000.nc test.nc
bsub -I -q q_x86_expr -n 2 ./test ../../../testdata/ocean/ocean_month.nc oceantest.nc
