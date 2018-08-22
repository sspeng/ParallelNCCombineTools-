/*************************************************************************
	> File Name: Argument.h
	> Author: ShiShupeng
	> Mail: shishupeng@mail.nsccwx.cn 
	> Created Time: 2018年08月01日 星期三 14时33分00秒
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

typedef struct fileinfo
{
	int filesnum;
    int ncfid;  /* ID of the input netCDF file */
    int ndims;  /* Number of dimensions */
    int nvars;  /* Number of variables */
    int ngatts;  /* Number of global attributes */
    int recdim;  /* ID of the record dimensions */
    char varname[MAX_NC_VARS][MAX_NC_NAME];  /* Names of the variables */
    nc_type datatype[MAX_NC_VARS]; /* Data types of the variables */
    int varndims[MAX_NC_VARS];  /* Number of dimensions for each variable */
    int vardim[MAX_NC_VARS][MAX_NC_DIMS];  /* Dimensions for each variable */
    int natts[MAX_NC_VARS];  /* Number of attributes for each variable */
	int varid[MAX_NC_VARS];
    unsigned char vardecomp[MAX_NC_VARS];  /* Is the variable decomposed */
    char dimname[MAX_NC_DIMS][MAX_NC_NAME];  /* Names of the dimensions */
    long dimsize[MAX_NC_DIMS];  /* Sizes of the dimensions (decomposed) */
    long dimfullsize[MAX_NC_DIMS];  /* Full sizes of the dimensions */
    long dimstart[MAX_NC_DIMS];  /* Start positions within full dimensions */
    long dimend[MAX_NC_DIMS];  /* End positions within full dimensions */
    unsigned char varmiss[MAX_NC_VARS];  /* Does variable have missing_value */
    unsigned char varmissval[MAX_NC_VARS][8];  /* missing_value per variable */
} FileInfo;

typedef struct AttInfo {
	char name[NC_MAX_NAME + 1];
	nc_type type;
	int len;
	char data[512];
	int datai[4];
	double datad[4];
	float dataf[4];

} Att;

typedef struct VarInfo {
	int id;
	char name[NC_MAX_NAME + 1];
	nc_type type;
	int ndims;
	int dimid[4];
	int dimfullsize[4];
	int start[4];
	int end[4];
	int attnum;
	Att att[4];
} Var;

typedef struct BufInfo {
	int id;
	char name[NC_MAX_NAME + 1];
	nc_type type;
	int ndims;
	int dimid[4];
	size_t varsize;
	int Xst;
	int Xed;
	int Yst;
	int Yed;
	int Zst;
	int Zed;
} Buf;

double *alldatad = NULL;
float *alldataf = NULL;
int *alldatai = NULL;
char *alldatac = NULL;
short *alldatas = NULL;

double *datad = NULL;
float *dataf = NULL;
int *datai = NULL;
char *datac = NULL;
short *datas = NULL;

void InitOutfile(FileInfo *ncoutfile, char *outncname)
{
	FileInfo *ncinfile;
	int nfiles2;
	int d, v, n;
	int dimid;
	int decomp[4];
	char attname[MAX_NC_NAME];
	
	ncinfile = (FileInfo*) malloc (sizeof(FileInfo));

	char filename[256];
	sprintf(filename, "%s.%04d", outncname, 0);

	ncinfile->ncfid = ncopen(filename, NC_NOWRITE);

	ncattget(ncinfile->ncfid, NC_GLOBAL, "NumFilesInSet", (void*) &nfiles2);

	ncoutfile->filesnum = nfiles2;

    /* Get some general information about the input netCDF file */
    if (ncinquire(ncinfile->ncfid,&(ncinfile->ndims),&(ncinfile->nvars),
                  &(ncinfile->ngatts),&(ncinfile->recdim))==(-1)) {
        fprintf(stderr,"Error: cannot read the file's metadata!\n");
        ncclose(ncinfile->ncfid); free(ncinfile); 
    }
	ncoutfile->ndims = ncinfile->ndims;
	ncoutfile->ngatts = ncinfile->ngatts;

    /* Get some information about the dimensions */
    for (d=0; d < ncinfile->ndims; d++) {
        if ((ncdiminq(ncinfile->ncfid,d,ncinfile->dimname[d],
                     &(ncinfile->dimsize[d])))==(-1)) {
            fprintf(stderr,"Error: cannot read dimension #%d's metadata!\n",d);
            ncclose(ncinfile->ncfid); free(ncinfile); 
        }
        ncinfile->dimfullsize[d]=ncinfile->dimsize[d];
        ncinfile->dimstart[d]=1; 
		ncinfile->dimend[d]=(-1);
    }

    /* Save some information for the output file */
	ncoutfile->nvars=ncinfile->nvars; 
	ncoutfile->recdim=ncinfile->recdim;

    /* Get some information about the variables */
    for (v=0; v < ncinfile->nvars; v++) {
        if ((ncvarinq(ncinfile->ncfid,v,ncinfile->varname[v],
                     &(ncinfile->datatype[v]),&(ncinfile->varndims[v]),
                     ncinfile->vardim[v],&(ncinfile->natts[v])))==(-1)) {
            fprintf(stderr,"Error: cannot read variable #%d's metadata!\n",v);
            ncclose(ncinfile->ncfid); 
			free(ncinfile); 
        }
		ncoutfile->datatype[v] = ncinfile->datatype[v];
		memcpy(ncoutfile->vardim[v], ncinfile->vardim[v], sizeof(int) * 4);
		ncoutfile->natts[v] = ncinfile->natts[v];

		dimid = -1;
		nc_inq_dimid(ncinfile->ncfid,ncinfile->varname[v], &dimid);
        /* If the variable is also a dimension then get decomposition info */
        if ((dimid)!=(-1)) {
			int jj = nc_get_att_int(ncinfile->ncfid,v,"domain_decomposition", (void *)decomp);
            if (jj == 0) {
                /* the dimension is decomposed */
                ncinfile->dimfullsize[dimid]=decomp[1]-decomp[0]+1;
                ncinfile->dimstart[dimid]=decomp[2]-(decomp[0]-1);
                ncinfile->dimend[dimid]=decomp[3]-(decomp[0]-1);
            }
            else {
				/* the dimension is NOT decomposed */
				ncinfile->dimfullsize[dimid]=ncinfile->dimsize[dimid];
				ncinfile->dimstart[dimid]=1; ncinfile->dimend[dimid]=(-1);
            }
        }
    }

    /* Get some additional information about the variables */
    for (v=0; v < ncinfile->nvars; v++) {

        /* start by assuming the variable has no decomposed dimension */
        ncinfile->vardecomp[v]=0;

        /* now, iterate over the variable's dimensions and mark the */
        /* variable as a decomposed variable if any dimension of */
        /* the variable is decomposed */
        for (d=0; d < ncinfile->varndims[v]; d++) {
            /* Does the variable have a decomposed dimension? */
            if (ncinfile->dimend[ncinfile->vardim[v][d]]!=(-1)) {
                ncinfile->vardecomp[v]=1; 
				break;
            }
        }

		//printf("varid%d: dim = %d\n", v, ncoutfile->varndims[v]);
        /* Save some information for the output file */
        /* This only needs to be done once per output file */
        ncoutfile->varndims[v]=ncinfile->varndims[v];
        for (d=0; d < ncinfile->varndims[v]; d++) {
            ncoutfile->vardim[v][d]=ncinfile->vardim[v][d];
		}

        ncoutfile->vardecomp[v]=ncinfile->vardecomp[v];
        strcpy(ncoutfile->varname[v],ncinfile->varname[v]);
        ncoutfile->varmiss[v]=0;
    }

    for (d=0; d < ncinfile->ndims; d++) {
        ncoutfile->dimfullsize[d]=ncinfile->dimfullsize[d];
	}
	
    /* If the output netCDF file was just created then define its structure */

    /* Define the dimensions */
    for (d=0; d < ncinfile->ndims; d++) {
		ncdimdef(ncoutfile->ncfid,ncinfile->dimname[d], ncinfile->dimfullsize[d]);
    }

    /* Define the variables and copy their attributes */
    for (v=0; v < ncinfile->nvars; v++) {
        ncoutfile->varid[v] = ncvardef(ncoutfile->ncfid,ncinfile->varname[v],ncinfile->datatype[v], ncinfile->varndims[v],ncinfile->vardim[v]);
        for (n=0; n < ncinfile->natts[v]; n++) {
            ncattname(ncinfile->ncfid,v,n,attname);
            if (!strcmp(attname,"domain_decomposition")) continue;
            else {
                if (ncattcopy(ncinfile->ncfid,v,attname,ncoutfile->ncfid,v)==(-1)) {
                    fprintf(stderr,"Error: cannot copy variable \"%s\"'s attributes!\n",
                       ncinfile->varname[v]);
                    free(ncinfile); 
                }
            }
        }
    }

    /* Copy the global attributes */
    for (n=0; n < ncinfile->ngatts; n++) {
        ncattname(ncinfile->ncfid,NC_GLOBAL,n,attname);
        if (!strcmp(attname,"NumFilesInSet")) {
			continue;
		}
        else if (!strcmp(attname,"filename")) {
            ncattput(ncoutfile->ncfid,NC_GLOBAL,attname,NC_CHAR,
                 strlen(outncname),(void *)outncname);
		}
        else {
            if (ncattcopy(ncinfile->ncfid,NC_GLOBAL,attname,ncoutfile->ncfid, NC_GLOBAL)==(-1)) {
                fprintf(stderr,"Error: cannot copy the file's global attributes!\n");
            }
        }
    }

    nc_enddef(ncoutfile->ncfid);
}

int GetNvars(char *name)
{
	char filename[256];
	sprintf(filename, "%s.%04d", name, 0);
    int ncid;
    int nvars, varid,  attnum, natts;
    XXnc(nc_open(filename, 0, &ncid));

    for (nvars = 0;; nvars++)
        if (nc_inq_varname(ncid, nvars, NULL))
            break;
	return nvars;
}

void GetHeadinfo(Var *ncvar, char *name, int nvars, Att *gatt, int *gattnum) 
{
	char filename[256];
	sprintf(filename, "%s.%04d", name, 0);
    int ncid;
    int varid,  attnum, natts;
    XXnc(nc_open(filename, 0, &ncid));

	for(varid = 0; varid < nvars; varid++) {
        char namevar[NC_MAX_NAME + 1];
        nc_type xxtype;
        int ndims, natts, attnum, dimid[4];
        XXnc(nc_inq_var(ncid, varid, namevar, &xxtype, &ndims, dimid, &natts));
		ncvar[varid].id = varid;
		ncvar[varid].type = xxtype;
		ncvar[varid].ndims = ndims;
		strcpy(ncvar[varid].name, namevar);
		memcpy(ncvar[varid].dimid, dimid, sizeof(int) * ndims);

        for (attnum = 0; ;attnum++)
        {
            char name[NC_MAX_NAME + 1];
            if (nc_inq_attname(ncid, varid, attnum, name))
                break;
            nc_type xtype;
            XXnc(nc_inq_atttype(ncid, varid, name, &xtype));
			strcpy(ncvar[varid].att[attnum].name, name);
			ncvar[varid].att[attnum].type = xtype;

            if (xtype == NC_CHAR)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                char v[len + 1];
                v[len] = '\0';
                XXnc(nc_get_att_text(ncid, varid, name, v));
				ncvar[varid].att[attnum].len = len;
				strcpy(ncvar[varid].att[attnum].data, v);
            }
            else if (xtype == NC_DOUBLE)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                double v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_double(ncid, varid, name, &v[kk]));
				}
            }
            else if (xtype == NC_INT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                int v[len];
				int kk;
				XXnc(nc_get_att_int(ncid, varid, name, v));
				for(kk = 0; kk < len; kk++) {
					//printf("INT  %d\n", v[kk]);
				}
            }
            else if (xtype == NC_FLOAT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                float v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_float(ncid, varid, name, &v[kk]));
				}
            }
            else
            {
                printf("null\n");
            }
        }
		ncvar[varid].attnum = attnum;

		/*------------------------
		 * global att
		 * ------------------------*/
        for (attnum = 0; ;attnum++)
        {
            char name[NC_MAX_NAME + 1];
            if (nc_inq_attname(ncid, NC_GLOBAL, attnum, name))
                break;
            nc_type xtype;
            XXnc(nc_inq_atttype(ncid, NC_GLOBAL, name, &xtype));
			strcpy(gatt[attnum].name, name);
			gatt[attnum].type = xtype;

            if (xtype == NC_CHAR)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, NC_GLOBAL, name, &len));
                char v[len + 1];
                v[len] = '\0';
                XXnc(nc_get_att_text(ncid, NC_GLOBAL, name, v));
				gatt[attnum].len = len;
				strcpy(gatt[attnum].data, v);
            }
            else if (xtype == NC_DOUBLE)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, NC_GLOBAL, name, &len));
                double v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_double(ncid, NC_GLOBAL, name, &v[kk]));
				}
				gatt[attnum].len = len;
				memcpy(gatt[attnum].datad, v, sizeof(double) * len);
            }
            else if (xtype == NC_INT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, NC_GLOBAL, name, &len));
                int v[len];
				int kk;
				XXnc(nc_get_att_int(ncid, NC_GLOBAL, name, v));
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_int(ncid, NC_GLOBAL, name, &v[kk]));
				}
				gatt[attnum].len = len;
				memcpy(gatt[attnum].datai, v, sizeof(int) * len);
            }
            else if (xtype == NC_FLOAT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, NC_GLOBAL, name, &len));
                float v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_float(ncid, NC_GLOBAL, name, &v[kk]));
				}
				gatt[attnum].len = len;
				memcpy(gatt[attnum].dataf, v, sizeof(float) * len);
            }
            else
            {
                printf("null\n");
            }
        }
		*gattnum = attnum;
	}
	nc_close(ncid);
}

void GetVarinfo_all(int nvars, char *name, size_t *dimlens, Var *ncvar) 
{
	int varid;
	int ncid;
	char filename[256];
	sprintf(filename, "%s.%04d", name, 0);
    XXnc(nc_open(filename, 0, &ncid));
	int i;
	for(varid = 0; varid < nvars; varid++) {
        char namevar[NC_MAX_NAME + 1];
        nc_type xxtype;
        int ndims, natts, attnum, dimid[4];
        XXnc(nc_inq_var(ncid, varid, namevar, &xxtype, &ndims, dimid, &natts));
		ncvar[varid].id = varid;
		ncvar[varid].type = xxtype;
		ncvar[varid].ndims = ndims;
		memcpy(ncvar[varid].dimid, dimid, sizeof(int) * ndims);
		strcpy(ncvar[varid].name, namevar);

        for (attnum = 0; ;attnum++)
        {
            char name[NC_MAX_NAME + 1];
            if (nc_inq_attname(ncid, varid, attnum, name))
                break;
            nc_type xtype;
            XXnc(nc_inq_atttype(ncid, varid, name, &xtype));
			memcpy(ncvar[varid].att[attnum].name, name, strlen(name));
			ncvar[varid].att[attnum].type = xtype;

            if (xtype == NC_CHAR)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                char v[len + 1];
                v[len] = '\0';
                XXnc(nc_get_att_text(ncid, varid, name, v));
				ncvar[varid].att[attnum].len = len;
				strcpy(ncvar[varid].att[attnum].data, v);
            }
            else if (xtype == NC_DOUBLE)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                double v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_double(ncid, varid, name, &v[kk]));
				}
            }
            else if (xtype == NC_INT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                int v[len];
				int kk;
				XXnc(nc_get_att_int(ncid, varid, name, v));
				for(kk = 0; kk < len; kk++) {
					//printf("INT  %d\n", v[kk]);
				}
            }
            else if (xtype == NC_FLOAT)
            {
                size_t len;
                XXnc(nc_inq_attlen(ncid, varid, name, &len));
                float v[len];
				int kk;
				for(kk = 0; kk < len; kk++) {
					XXnc(nc_get_att_float(ncid, varid, name, &v[kk]));
				}
            }
            else
            {
                printf("null\n");
            }
        }
		ncvar[varid].attnum = attnum;
	}
}

void GetBufinfo(int ncid, int varid, size_t *dimlens, Buf *ncvar) 
{
    char namevar[NC_MAX_NAME + 1];
    nc_type xxtype;
    int ndims, natts, attnum, dimid[4];
    XXnc(nc_inq_var(ncid, varid, namevar, &xxtype, &ndims, dimid, &natts));
    ncvar[varid].id = varid;
    ncvar[varid].type = xxtype;
    ncvar[varid].ndims = ndims;
    memcpy(ncvar[varid].dimid, dimid, sizeof(int) * ndims);
    strcpy(ncvar[varid].name, namevar);

	int k;
	size_t varsize = 1;
	for(k = 0; k < ncvar[varid].ndims; k++) {
		varsize *= dimlens[ncvar[varid].dimid[k]];
	}
}

void GetVarinfo(int ncid, int varid, size_t *dimlens, Var *ncvar) 
{
    char namevar[NC_MAX_NAME + 1];
    nc_type xxtype;
    int ndims, natts, attnum, dimid[4];
    XXnc(nc_inq_var(ncid, varid, namevar, &xxtype, &ndims, dimid, &natts));
    ncvar->id = varid;
    ncvar->type = xxtype;
    ncvar->ndims = ndims;
    memcpy(ncvar->dimid, dimid, sizeof(int) * ndims);
    strcpy(ncvar->name, namevar);
}

void WriteData(int ncid, int varid, size_t varsize, Var ncvar, size_t *dimlens, size_t *dimlensall)
{
	size_t i, j, k, l, m, n, p, q;

	if(ncvar.ndims == 1) {
		double t1, t2, t3;
		char dimnames[256];
		int decomp[4];
		int dimvarid;
		size_t start, end, dimsize, dimfullsize;
		
        XXnc(nc_inq_dimname(ncid, ncvar.dimid[0], dimnames));
        XXnc(nc_inq_varid(ncid, dimnames, &dimvarid));
		int jj = nc_get_att(ncid, dimvarid, "domain_decomposition", (void *)decomp);

		if(jj == 0) {
			dimfullsize = decomp[1] - decomp[0] + 1;
			start = decomp[2] - decomp[0] + 1;
			end = decomp[3] - decomp[0] + 1;
			dimsize = end - start + 1;
		}
		else {
			start = 1;
			end = dimlensall[ncvar.dimid[0]];
			dimsize = end - start + 1;
		}
		
		switch (ncvar.type) { 
			case NC_BYTE:
			case NC_CHAR:
                #pragma omp parallel for 
				for(i = 0; i < dimsize; i++) {
					alldatac[i + start - 1] = datac[i]; 
				}
				break;
			case NC_SHORT:
                #pragma omp parallel for 
				for(i = 0; i < dimsize; i++) {
					alldatas[i + start - 1] = datas[i]; 
				}
				break;
			case NC_INT:
                #pragma omp parallel for 
				for(i = 0; i < dimsize; i++) {
					alldatai[i + start - 1] = datai[i]; 
				}
				break;
			case NC_FLOAT:
                #pragma omp parallel for 
				for(i = 0; i < dimsize; i++) {
					alldataf[i + start - 1] = dataf[i]; 
				}
				break;
			case NC_DOUBLE:
                #pragma omp parallel for 
				for(i = 0; i < dimsize; i++) {
					alldatad[i + start - 1] = datad[i]; 
				}
				break;
		}
	} else if(ncvar.ndims == 2) {
		double t1, t2, t3;
		char dimnames[2][256];
		int dimvarid[2];
		size_t start[2], end[2], dimsize[2], dimfullsize[2];
		int decomp[2][4];
		
        XXnc(nc_inq_dimname(ncid, ncvar.dimid[0], dimnames[0]));
        XXnc(nc_inq_varid(ncid, dimnames[0], &dimvarid[0]));
		int ii = nc_get_att(ncid, dimvarid[0], "domain_decomposition", (void *)decomp[0]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[1], dimnames[1]));
        XXnc(nc_inq_varid(ncid, dimnames[1], &dimvarid[1]));
		int jj = nc_get_att(ncid, dimvarid[1], "domain_decomposition", (void *)decomp[1]);

		if(ii == 0) {
			dimfullsize[0] = decomp[0][1] - decomp[0][0] + 1;
			start[0] = decomp[0][2] - decomp[0][0] + 1;
			end[0] = decomp[0][3] - decomp[0][0] + 1;
			dimsize[0] = end[0] - start[0] + 1;
		}
		else {
			start[0] = 1;
			end[0] = dimlensall[ncvar.dimid[0]];
			dimsize[0] = end[0] - start[0] + 1;
			dimfullsize[0] = dimsize[0];
		}


		if(jj == 0) {
			dimfullsize[1] = decomp[1][1] - decomp[1][0] + 1;
			start[1] = decomp[1][2] - decomp[1][0] + 1;
			end[1] = decomp[1][3] - decomp[1][0] + 1;
			dimsize[1] = end[1] - start[1] + 1;
		}
		else {
			start[1] = 1;
			end[1] = dimlensall[ncvar.dimid[1]];
			dimsize[1] = end[1] - start[1] + 1;
			dimfullsize[1] = dimsize[1];
		}

		switch (ncvar.type) { 
			case NC_BYTE:
			case NC_CHAR:
				#pragma omp parallel for private(i, j, k, l, m)
				for(j = 0; j < dimsize[0]; j++) {
					for(i = 0; i < dimsize[1]; i++) {
						l = j + start[0] - 1;
						k = i + start[1] - 1;
						m = l * dimfullsize[1] + k;
						alldatac[m] = datac[j * dimsize[1] + i]; 
					}
				}
				break;
			case NC_SHORT:
				#pragma omp parallel for private(i, j, k, l, m)
				for(j = 0; j < dimsize[0]; j++) {
					for(i = 0; i < dimsize[1]; i++) {
						l = j + start[0] - 1;
						k = i + start[1] - 1;
						m = l * dimfullsize[1] + k;
						alldatas[m] = datas[j * dimsize[1] + i]; 
					}
				}
				break;
			case NC_INT:
				#pragma omp parallel for private(i, j, k, l, m)
				for(j = 0; j < dimsize[0]; j++) {
					for(i = 0; i < dimsize[1]; i++) {
						l = j + start[0] - 1;
						k = i + start[1] - 1;
						m = l * dimfullsize[1] + k;
						alldatai[m] = datai[j * dimsize[1] + i]; 
					}
				}
				break;
			case NC_FLOAT:
				#pragma omp parallel for private(i, j, k, l, m)
				for(j = 0; j < dimsize[0]; j++) {
					for(i = 0; i < dimsize[1]; i++) {
						l = j + start[0] - 1;
						k = i + start[1] - 1;
						m = l * dimfullsize[1] + k;
						alldataf[m] = dataf[j * dimsize[1] + i]; 
					}
				}
				break;
			case NC_DOUBLE:
				#pragma omp parallel for private(i, j, k, l, m)
				for(j = 0; j < dimsize[0]; j++) {
					for(i = 0; i < dimsize[1]; i++) {
						l = j + start[0] - 1;
						k = i + start[1] - 1;
						m = l * dimfullsize[1] + k;
						alldatad[m] = datad[j * dimsize[1] + i]; 
					}
				}
				break;
		}
	} else if (ncvar.ndims == 3) {
		char dimnames[3][256];
		int dimvarid[3];
		size_t start[3], end[3], dimsize[3], dimfullsize[3];
		int decomp[3][4];
		
        XXnc(nc_inq_dimname(ncid, ncvar.dimid[0], dimnames[0]));
        XXnc(nc_inq_varid(ncid, dimnames[0], &dimvarid[0]));
		int ii = nc_get_att(ncid, dimvarid[0], "domain_decomposition", (void *)decomp[0]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[1], dimnames[1]));
        XXnc(nc_inq_varid(ncid, dimnames[1], &dimvarid[1]));
		int jj = nc_get_att(ncid, dimvarid[1], "domain_decomposition", (void *)decomp[1]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[2], dimnames[2]));
        XXnc(nc_inq_varid(ncid, dimnames[2], &dimvarid[2]));
		int kk = nc_get_att(ncid, dimvarid[2], "domain_decomposition", (void *)decomp[2]);

		if(ii == 0) {
			dimfullsize[0] = decomp[0][1] - decomp[0][0] + 1;
			start[0] = decomp[0][2] - decomp[0][0] + 1;
			end[0] = decomp[0][3] - decomp[0][0] + 1;
			dimsize[0] = end[0] - start[0] + 1;
		}
		else {
			start[0] = 1;
			end[0] = dimlensall[ncvar.dimid[0]];
			dimsize[0] = end[0] - start[0] + 1;
			dimfullsize[0] = dimsize[0];
		}

		if(jj == 0) {
			dimfullsize[1] = decomp[1][1] - decomp[1][0] + 1;
			start[1] = decomp[1][2] - decomp[1][0] + 1;
			end[1] = decomp[1][3] - decomp[1][0] + 1;
			dimsize[1] = end[1] - start[1] + 1;
		}
		else {
			start[1] = 1;
			end[1] = dimlensall[ncvar.dimid[1]];
			dimsize[1] = end[1] - start[1] + 1;
			dimfullsize[1] = dimsize[1];
		}

		if(kk == 0) {
			dimfullsize[2] = decomp[2][1] - decomp[2][0] + 1;
			start[2] = decomp[2][2] - decomp[2][0] + 1;
			end[2] = decomp[2][3] - decomp[2][0] + 1;
			dimsize[2] = end[2] - start[2] + 1;
		}
		else {
			start[2] = 1;
			end[2] = dimlensall[ncvar.dimid[2]];
			dimsize[2] = end[2] - start[2] + 1;
			dimfullsize[2] = dimsize[2];
		}

		size_t jks = dimsize[2] * dimsize[1];
		size_t jksall = dimfullsize[2] * dimfullsize[1];

		size_t ks = dimsize[2];
		size_t ksall = dimfullsize[2];

		switch (ncvar.type) { 
			case NC_BYTE:
			case NC_CHAR:
				#pragma omp parallel for private(i, j, k, l, m, n)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							m = k + start[0] - 1;
							n = j + start[1] - 1;
							l = i + start[2] - 1;

							alldatac[m * jksall + n * ksall + l] = datac[k * jks + j * ks + i]; 
						}
					}
				}
				break;
			case NC_SHORT:
				#pragma omp parallel for private(i, j, k, l, m, n)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							m = k + start[0] - 1;
							n = j + start[1] - 1;
							l = i + start[2] - 1;

							alldatas[m * jksall + n * ksall + l] = datas[k * jks + j * ks + i]; 
						}
					}
				}
				break;
			case NC_INT:
				#pragma omp parallel for private(i, j, k, l, m, n)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							m = k + start[0] - 1;
							n = j + start[1] - 1;
							l = i + start[2] - 1;

							alldatai[m * jksall + n * ksall + l] = datai[k * jks + j * ks + i]; 
						}
					}
				}
				break;
			case NC_FLOAT:
				#pragma omp parallel for private(i, j, k, l, m, n)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							m = k + start[0] - 1;
							n = j + start[1] - 1;
							l = i + start[2] - 1;

							alldataf[m * jksall + n * ksall + l] = dataf[k * jks + j * ks + i]; 
						}
					}
				}
				break;
			case NC_DOUBLE:
				#pragma omp parallel for private(i, j, k, l, m, n)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							m = k + start[0] - 1;
							n = j + start[1] - 1;
							l = i + start[2] - 1;

							alldatad[m * jksall + n * ksall + l] = datad[k * jks + j * ks + i]; 
						}
					}
				}
				break;
		}
	} else if (ncvar.ndims == 4){
		char dimnames[4][256];
		int dimvarid[4];
		size_t start[4], end[4], dimsize[4], dimfullsize[4];
		int decomp[4][4];
		
        XXnc(nc_inq_dimname(ncid, ncvar.dimid[0], dimnames[0]));
        XXnc(nc_inq_varid(ncid, dimnames[0], &dimvarid[0]));
		int ii = nc_get_att(ncid, dimvarid[0], "domain_decomposition", (void *)decomp[0]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[1], dimnames[1]));
        XXnc(nc_inq_varid(ncid, dimnames[1], &dimvarid[1]));
		int jj = nc_get_att(ncid, dimvarid[1], "domain_decomposition", (void *)decomp[1]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[2], dimnames[2]));
        XXnc(nc_inq_varid(ncid, dimnames[2], &dimvarid[2]));
		int kk = nc_get_att(ncid, dimvarid[2], "domain_decomposition", (void *)decomp[2]);

        XXnc(nc_inq_dimname(ncid, ncvar.dimid[3], dimnames[3]));
        XXnc(nc_inq_varid(ncid, dimnames[3], &dimvarid[3]));
		int ll = nc_get_att(ncid, dimvarid[3], "domain_decomposition", (void *)decomp[3]);

		if(ii == 0) {
			dimfullsize[0] = decomp[0][1] - decomp[0][0] + 1;
			start[0] = decomp[0][2] - decomp[0][0] + 1;
			end[0] = decomp[0][3] - decomp[0][0] + 1;
			dimsize[0] = end[0] - start[0] + 1;
		}
		else {
			start[0] = 1;
			end[0] = dimlensall[ncvar.dimid[0]];
			dimsize[0] = end[0] - start[0] + 1;
			dimfullsize[0] = dimsize[0];
		}

		if(jj == 0) {
			dimfullsize[1] = decomp[1][1] - decomp[1][0] + 1;
			start[1] = decomp[1][2] - decomp[1][0] + 1;
			end[1] = decomp[1][3] - decomp[1][0] + 1;
			dimsize[1] = end[1] - start[1] + 1;
		}
		else {
			start[1] = 1;
			end[1] = dimlensall[ncvar.dimid[1]];
			dimsize[1] = end[1] - start[1] + 1;
			dimfullsize[1] = dimsize[1];
		}

		if(kk == 0) {
			dimfullsize[2] = decomp[2][1] - decomp[2][0] + 1;
			start[2] = decomp[2][2] - decomp[2][0] + 1;
			end[2] = decomp[2][3] - decomp[2][0] + 1;
			dimsize[2] = end[2] - start[2] + 1;
		}
		else {
			start[2] = 1;
			end[2] = dimlensall[ncvar.dimid[2]];
			dimsize[2] = end[2] - start[2] + 1;
			dimfullsize[2] = dimsize[2];
		}

		if(ll == 0) {
			dimfullsize[3] = decomp[3][1] - decomp[3][0] + 1;
			start[3] = decomp[3][2] - decomp[3][0] + 1;
			end[3] = decomp[3][3] - decomp[3][0] + 1;
			dimsize[3] = end[3] - start[3] + 1;
		}
		else {
			start[3] = 1;
			end[3] = dimlensall[ncvar.dimid[3]];
			dimsize[3] = end[3] - start[3] + 1;
			dimfullsize[3] = dimsize[3];
		}

		size_t jils = dimsize[3] * dimsize[2] * dimsize[1];
		size_t jilsall = dimfullsize[3] * dimfullsize[2] * dimfullsize[1];

		size_t ils = dimsize[3] * dimsize[2];
		size_t ilsall = dimfullsize[3] * dimfullsize[2];

		size_t ls = dimsize[3];
		size_t lsall = dimfullsize[3];

		switch (ncvar.type) { 
			case NC_BYTE:
			case NC_CHAR:
				#pragma omp parallel for private(k, l, m, n, p, q)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							for(l = 0; l < dimsize[3]; l++) {
								m = k + start[0] - 1;
								n = j + start[1] - 1;
								p = i + start[2] - 1;
								q = l + start[3] - 1;
								alldatac[m * jilsall + n * ilsall + p * lsall + q] = datac[k * jils + j * ils + i * ls + l]; 
							}
						}
					}
				}
				break;
			case NC_SHORT:
				#pragma omp parallel for private(k, l, m, n, p, q)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							for(l = 0; l < dimsize[3]; l++) {
								m = k + start[0] - 1;
								n = j + start[1] - 1;
								p = i + start[2] - 1;
								q = l + start[3] - 1;
								alldatas[m * jilsall + n * ilsall + p * lsall + q] = datas[k * jils + j * ils + i * ls + l]; 
							}
						}
					}
				}
				break;
			case NC_INT:
				#pragma omp parallel for private(k, l, m, n, p, q)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							for(l = 0; l < dimsize[3]; l++) {
								m = k + start[0] - 1;
								n = j + start[1] - 1;
								p = i + start[2] - 1;
								q = l + start[3] - 1;
								alldatai[m * jilsall + n * ilsall + p * lsall + q] = datai[k * jils + j * ils + i * ls + l]; 
							}
						}
					}
				}
				break;
			case NC_FLOAT:
				#pragma omp parallel for private(k, l, m, n, p, q)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							for(l = 0; l < dimsize[3]; l++) {
								m = k + start[0] - 1;
								n = j + start[1] - 1;
								p = i + start[2] - 1;
								q = l + start[3] - 1;
								alldataf[m * jilsall + n * ilsall + p * lsall + q] = dataf[k * jils + j * ils + i * ls + l]; 
							}
						}
					}
				}
				break;
			case NC_DOUBLE:
				#pragma omp parallel for private(k, l, m, n, p, q)
				for(k = 0; k < dimsize[0]; k++) {
					for(j = 0; j < dimsize[1]; j++) {
						for(i = 0; i < dimsize[2]; i++) {
							for(l = 0; l < dimsize[3]; l++) {
								m = k + start[0] - 1;
								n = j + start[1] - 1;
								p = i + start[2] - 1;
								q = l + start[3] - 1;
								alldatad[m * jilsall + n * ilsall + p * lsall + q] = datad[k * jils + j * ils + i * ls + l]; 
							}
						}
					}
				}
				break;
		}
	} else {
		printf("ERROR: dims > 4\n");
		exit(-1);
	}
}
