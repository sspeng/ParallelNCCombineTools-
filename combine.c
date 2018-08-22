/*************************************************************************
	> File Name: para_com.c
	> Author: ShiShupeng
	> Mail: shishupeng@mail.nsccwx.cn 
	> Created Time: 2018年07月28日 星期六 11时06分56秒
	> Program Infomation: 
	> 
 ************************************************************************/
#include "Argument.h"

int main(int argc, char **argv)
{
	assert(argc == 3);

	int i, j, k, l, num, varid;
	int rank, size;
	int start, end, sec, res, mysec;
	double *ldoflen;
	double t1, t2;
	ldoflen = (double*) malloc(num * sizeof(double));

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	t1 = MPI_Wtime();

	int ncidall;
	nc_create_par(argv[2], NC_CLOBBER | NC_PNETCDF | NC_64BIT_OFFSET, MPI_COMM_WORLD, MPI_INFO_NULL, &ncidall);

	nc_set_fill(ncidall, NC_NOFILL, NULL);

	FileInfo *ncoutfile = (FileInfo*) malloc (sizeof(FileInfo));

	ncoutfile->ncfid = ncidall;
	
	InitOutfile(ncoutfile, argv[1]);	

	if(rank == 0) { 
		printf("Total Num: %d    nvars:  %d\n", ncoutfile->filesnum, ncoutfile->nvars);
	}

#ifdef POP_Balance
	int st, ed;
	if(rank == 0 && size != 14) {
		printf("Error: POP_Balance needs 14 MPI processes\n");
		MPI_Finalize();
		return 0;
	}
	if(rank == 0) {
		st = 0;
		ed = 5;
		printf("POP_Balance!\n");
	}
	else if(rank == 1) {
		st = 5;
		ed = 9;
	}
	else if(rank == 2) {
		st = 9;
		ed = 13;
	}
	else if(rank == 3) {
		st = 13;
		ed = 17;
	}
	else {
		st = 13 + rank;
		ed = 13 + rank + 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(varid = st; varid < ed; varid++) {
#else
	MPI_Barrier(MPI_COMM_WORLD);
	for(varid = rank; varid < ncoutfile->nvars; varid+=size) {
#endif
		size_t k, varsize = 1;
		int tag = 0;

		for(k = 0; k < ncoutfile->varndims[varid]; k++) {
            varsize*=ncoutfile->dimfullsize[ncoutfile->vardim[varid][k]];
		}

		nc_var_par_access(ncidall, varid, 0);

		switch (ncoutfile->datatype[varid]) { 
			case NC_BYTE:
			case NC_CHAR:
				alldatac = (char*) malloc (sizeof(char) * varsize);
				break;
			case NC_SHORT:
				alldatas = (short*) malloc (sizeof(short) * varsize);
				break;
			case NC_INT:
				alldatai = (int*) malloc (sizeof(int) * varsize);
				break;
			case NC_FLOAT:
				alldataf = (float*) malloc (sizeof(float) * varsize);
				break;
			case NC_DOUBLE:
				alldatad = (double*) malloc (sizeof(double) * varsize);
				break;
		}
		
		for(i = 0; i < ncoutfile->filesnum; i++) {
			int ncid;
			char filename[256];
			sprintf(filename, "%s.%04d", argv[1], i);

			nc_open(filename, 0, &ncid);

			int ndims;
			XXnc(nc_inq_ndims(ncid, &ndims));

			size_t dimlens[ndims];

			for (j = 0; j < ndims; j++) {
				XXnc(nc_inq_dimlen(ncid, j, &dimlens[j]));
			}

			Var ncvar;
			GetVarinfo(ncid, varid, dimlens, &ncvar); 

			int k;
			size_t varsize = 1;

			for(k = 0; k < ncvar.ndims; k++) {
				varsize *= dimlens[ncvar.dimid[k]];
			}

			switch (ncvar.type) { 
				case NC_BYTE:
				case NC_CHAR:
					datac = (char*) malloc (sizeof(char) * varsize);
					nc_get_var (ncid, varid, datac);
					WriteData (ncid, varid, varsize, ncvar, dimlens, ncoutfile->dimfullsize);
					free(datac);
					break;
				case NC_SHORT:
					datas = (short*) malloc (sizeof(short) * varsize);
					nc_get_var (ncid, varid, datas);
					WriteData (ncid, varid, varsize, ncvar, dimlens, ncoutfile->dimfullsize);
					free(datas);
					break;
				case NC_INT:
					datai = (int*) malloc (sizeof(int) * varsize);
					nc_get_var (ncid, varid, datai);
					WriteData (ncid, varid, varsize, ncvar, dimlens, ncoutfile->dimfullsize);
					free(datai);
					break;
				case NC_FLOAT:
					dataf = (float*) malloc (sizeof(float) * varsize);
					nc_get_var (ncid, varid, dataf);
					WriteData (ncid, varid, varsize, ncvar, dimlens, ncoutfile->dimfullsize);
					free(dataf);
					break;
				case NC_DOUBLE:
					datad = (double*) malloc (sizeof(double) * varsize);
					nc_get_var (ncid, varid, datad);
					WriteData (ncid, varid, varsize, ncvar, dimlens, ncoutfile->dimfullsize);
					free(datad);
					break;
			}
			
			nc_close(ncid);

			if ((ncoutfile->filesnum > 100) && ((i + 1) % (ncoutfile->filesnum / 10)) == 0) {
				tag += 1;
				printf("Rank%d:\t varid: %d\t %d%%\n", rank, varid, tag * 10);
			} 
		}

        switch (ncoutfile->datatype[varid]) { 
            case NC_BYTE:
            case NC_CHAR:
				nc_sync(ncidall);
				XXnc(nc_put_var_text(ncoutfile->ncfid, ncoutfile->varid[varid], alldatac));
				nc_sync(ncidall);
				free(alldatac);
                break;
		    case NC_SHORT:
				nc_sync(ncidall);
				XXnc(nc_put_var_short(ncoutfile->ncfid, ncoutfile->varid[varid], alldatas));
				nc_sync(ncidall);
				free(alldatas);
                break;
		    case NC_INT:
				nc_sync(ncidall);
				XXnc(nc_put_var_int(ncoutfile->ncfid, ncoutfile->varid[varid], alldatai));
				nc_sync(ncidall);
				free(alldatai);
                break;
		    case NC_FLOAT:
				nc_sync(ncidall);
				XXnc(nc_put_var_float(ncoutfile->ncfid, ncoutfile->varid[varid], alldataf));
				nc_sync(ncidall);
				free(alldataf);
                break;
		    case NC_DOUBLE:
				nc_sync(ncidall);
				XXnc(nc_put_var_double(ncoutfile->ncfid, ncoutfile->varid[varid], alldatad));
				nc_sync(ncidall);
				free(alldatad);
                break;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	nc_close(ncidall);

	t2 = MPI_Wtime();

	if(rank == size - 1)
	printf("Total Time: %lfmin\n", (t2 - t1)/60.);
	MPI_Finalize();

	return 0;
}
