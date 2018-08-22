/*************************************************************************
	> File Name: resultcheck.c
	> Author: ShiShupeng
	> Mail: shishupeng@mail.nsccwx.cn 
	> Created Time: 2018年08月14日 星期二 16时42分16秒
	> Program Infomation: 
	> 
 ************************************************************************/
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>

#define XXnc(cmd)                     \
{                                     \
    int status = (cmd);               \
    if (status)                       \
    {                                 \
        fprintf(stderr,               \
                "%s:%i: error: %s\n", \
                __FILE__,             \
                __LINE__,             \
                nc_strerror(status)); \
    }                                 \
}

int main(int argc, char **argv)
{
	assert(argc == 3);

	int i, j, k, l, num, varid;
	int rank, size;
	int start, end, sec, res, mysec;
	double t1, t2;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	t1 = MPI_Wtime();

	int ncid1, ncid2;

	nc_open(argv[1], NC_NOWRITE, &ncid1);
	nc_open(argv[2], NC_NOWRITE, &ncid2);

	int ndimsall;
	XXnc(nc_inq_ndims(ncid1, &ndimsall));

	size_t dimlens[ndimsall];

	for (j = 0; j < ndimsall; j++) {
		XXnc(nc_inq_dimlen(ncid1, j, &dimlens[j]));
	}

	int nvars;
    for (nvars = 0;; nvars++)
        if (nc_inq_varname(ncid1, nvars, NULL))
            break;

	for(varid = rank; varid < nvars; varid+=size) {

		int ndims;
		nc_type xxtype;
		int natts, attnum, dimid[4];
		char namevar[NC_MAX_NAME + 1];

		XXnc(nc_inq_var(ncid1, varid, namevar, &xxtype, &ndims, dimid, &natts));

		size_t k, varsize = 1;

		for(k = 0; k < ndims; k++) {
			varsize *= dimlens[dimid[k]];
		}

		double *vardata = (double*) malloc (sizeof(double) * varsize);
		double *vardata2 = (double*) malloc (sizeof(double) * varsize);

		int cvarid;
        XXnc(nc_inq_varid(ncid2, namevar, &cvarid));

		nc_get_var_double(ncid1, varid, vardata);
		nc_get_var_double(ncid2, cvarid, vardata2);

		for(i = 0; i < varsize; i++) {
			if(fabs(vardata[i] - vardata2[i]) > 1e-12) {
				printf("Error: varid=%d, varsize = %ld\n", varid, i);
			}
		}

		free(vardata);
		free(vardata2);

		printf("VARID%3d:\t %20s\t Result Correct!\n", varid, namevar);
	}

	nc_close(ncid1);
	nc_close(ncid2);

	t2 = MPI_Wtime();

	if(rank == size - 1)
	printf("Total Time: %lfmin\n", (t2 - t1)/60.);
	MPI_Finalize();

	return 0;
}

