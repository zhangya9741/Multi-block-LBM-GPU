
#include "parameters.h"

#include "geo/init_geo.h"

#include "lbm/assit_func.h"
#include "lbm/init_block.h"
#include "lbm/force_block.h"
#include "lbm/error_L2.h"

#include "multi/init_multi.h"
#include "multi/FH/assit_multi.h"
#include "multi/FH/evolution_block.h"

#include "output/output.h"

#include "evolution.h"

int main()
{
	geo_init();

	Lbui N_varible;
	D2Q9 Variables_h[number_area], Variables_d[number_area];
	/*To allocate lbm variables*/
	for (Lbi a=0; a<number_area; a++)
	{
		N_varible = Areas[a].nx*Areas[a].ny;
		Variables_h[a].solid = (Lbi *)malloc(N_varible*sizeof(Lbi));
		cudaMalloc((void **)&Variables_d[a].solid, N_varible*sizeof(Lbi));
		Lbd* *Vh = &Variables_h[a].x;
		Lbd* *Vd = &Variables_d[a].x;
		for (Lbi b=0; b<1+3*DMS+2*Q; b++)
		{
			Vh[b] = (Lbd *)malloc(N_varible*sizeof(Lbd));
			cudaMalloc((void **)&Vd[b], N_varible*sizeof(Lbd)); 
		}
	}

	/*To allocate multi-block variables*/
	/*Temperal interpolation*/
	Ff_Store Line_T3_C[number_line>0 ? number_line : 1];
	for (Lbi a=0; a<number_line; a++)
	{ 
		Line_T3_C[a].marker[0] = 2;
		Line_T3_C[a].marker[1] = 1;
		Line_T3_C[a].marker[2] = 0;

		Line_T3_C[a].N = Line_Inf_C[a].Nf;

		F_k3 *ft = &Line_T3_C[a].f_t1;
		for (Lbi i=0; i<3; i++) //3t
		{
			Lbd **ftk = &ft[i].F[0];
			for (Lbi j=0; j<Q; j++) //Q
				cudaMalloc((void **)&ftk[j], Line_T3_C[a].N*sizeof(Lbd));
		}
	}

	cusparseHandle_t handle;
	cusparseCreate(&handle);

	/*Spacial interpolation*/
	FM_spline3 Line_Spline3_C[number_line>0 ? number_line : 1];
	for (Lbi a=0; a<number_line; a++)
	{
		Line_Spline3_C[a].N = Line_Inf_C[a].Nc;

		cudaMalloc((void **)&Line_Spline3_C[a].l, Line_Spline3_C[a].N*sizeof(Lbd));
		cudaMalloc((void **)&Line_Spline3_C[a].m, Line_Spline3_C[a].N*sizeof(Lbd));
		cudaMalloc((void **)&Line_Spline3_C[a].u, Line_Spline3_C[a].N*sizeof(Lbd));
		cudaMalloc((void **)&Line_Spline3_C[a].Varibs_M, 12*Line_Spline3_C[a].N*sizeof(Lbd));

		Lbd* *fM = &Line_Spline3_C[a].rho;
		for (Lbi i=0; i<12; i++)
			cudaMalloc((void **)&fM[i], Line_Spline3_C[a].N*sizeof(Lbd));
	}

	/*Get the block size of GPU kernels of LBM evolution*/
	dim3 block(_BLOCKX, _BLOCKY, 1);
	dim3 *grid = (dim3 *)malloc(number_area*sizeof(dim3));
	for (Lbi i=0; i<number_area; i++)
	{
		dim3 grid_tmp(int (ceil(Areas[i].nx/double (_BLOCKX))), int (ceil(Areas[i].ny/double (_BLOCKY))), 1);
		grid[i] = grid_tmp;
	} 

	/*To store the fluid velocity at the previous time step*/
	Lbd *upsx[number_area], *upsy[number_area];

	Lbd *h_uvp2_b[number_area];
	Lbd *d_uvp2_b[number_area];
	Lbd *h_uv2_b[number_area];
	Lbd *d_uv2_b[number_area];

	for (Lbi a=0; a<number_area; a++)
	{
		N_varible = Areas[a].nx*Areas[a].ny;
		cudaMalloc((void **)&upsx[a], N_varible*sizeof(Lbd));
		cudaMalloc((void **)&upsy[a], N_varible*sizeof(Lbd));

		Lbi Bx = int (ceil(Areas[a].nx/double (_BLOCKX)));
		Lbi By = int (ceil(Areas[a].ny/double (_BLOCKY)));
		Lbui size_B = Bx*By*sizeof(Lbd);

		h_uvp2_b[a] = (Lbd *)malloc(size_B);
		d_uvp2_b[a] = NULL;
		cudaMalloc((void **)&d_uvp2_b[a], size_B);

		h_uv2_b[a] = (Lbd *)malloc(size_B);
		d_uv2_b[a] = NULL;
		cudaMalloc((void **)&d_uv2_b[a], size_B);
	}

	/*To initialize the computational domain*/
	for (Lbi a=0; a<number_area; a++)
		init_General<<<grid[a], block>>>(Areas[a], Variables_d[a], upsx[a], upsy[a]);

	/*Output lbm variables*/
	output(0, Areas, Variables_h, Variables_d);
	
	/************************************************************************/

	/*To allocate variables for calculating the darg force with the momentum-exchange (ME) method*/
	Lbd Force[2];
	//Return the total number of lattices taking part in the ME
	Lbi N_drag = sum_N_drag(block, grid[number_area-1], Areas[number_area-1], Variables_d[number_area-1]);
	Int3 *d_Point_drag;
	cudaMalloc((void **)&d_Point_drag, N_drag*sizeof(Int3));
	//Store the location (i, j ~ LX, LY) and the discrete velocity (k~Q) to ease the following drag force calculation
	fdrag_store(Areas[number_area-1], Variables_h[number_area-1], d_Point_drag, Variables_d[number_area-1], N_drag);
	
	/************************************************************************/

	/*To initialize the variables of multi-block*/
	for (Lbi a=0; a<number_line; a++)
	{
		init_Fstore<<< int (ceil(Line_T3_C[a].N/double (_BLOCK_1D))), _BLOCK_1D>>>(Areas[Line_Inf_C[a].areaf_n], Variables_d[Line_Inf_C[a].areaf_n], Line_Inf_C[a], Line_T3_C[a]);
		init_tri_A<<< int (ceil(Line_Spline3_C[a].N/double (_BLOCK_1D))), _BLOCK_1D>>>(Areas[Line_Inf_C[a].areac_n], Variables_d[Line_Inf_C[a].areac_n], Line_Inf_C[a], Line_Spline3_C[a]);
	}

	clock_t start, stop;
	Lbd duration;
	for (Lbi k=0; k<_LOOPS; k++)
	{
		duration = 0.0;
		start = clock();
		for (Lbi n=0; n<_T; n++)
		{
			/************************************************************************/
			evolution(k*_T+n, ASP, LS3CFFCP, LTCFFCP, LCCP, 0, block, grid, Areas, Variables_h, 
				Variables_d, Line_Inf_C, Line_Inf_F, Line_Inf_CC, Line_T3_C, Line_Spline3_C, handle);
			/************************************************************************/

			error_L2(k*_T+n, number_area, block, grid, Areas, Variables_d, upsx, upsy, h_uvp2_b, d_uvp2_b, h_uv2_b, d_uv2_b);

			getDrag_block(k*_T+n, d_Point_drag, Areas[number_area-1], Variables_d[number_area-1], Force, N_drag);
		}

		stop = clock();
		duration = (Lbd)(stop-start)/CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);

		FILE *fp_t;
		fp_t=fopen("time_t.dat","a");
		fprintf(fp_t,"%d %g\n", (k+1)*_T, duration);
		fclose(fp_t);

		printf("%d\n", (k+1)*_T);
		output((k+1)*_T, Areas, Variables_h, Variables_d);
	}

	for (Lbi a=0; a<number_line; a++)
	{
		cudaFree(Line_Spline3_C[a].l);
		cudaFree(Line_Spline3_C[a].m);
		cudaFree(Line_Spline3_C[a].u);
		cudaFree(Line_Spline3_C[a].Varibs_M);
		Lbd* *fM = &Line_Spline3_C[a].rho;
		for (Lbi i=0; i<12; i++)
			cudaFree(fM[i]);
	}
	cusparseDestroy(handle);

	for (Lbi a=0; a<number_line; a++)
	{ 
		F_k3 *ft = &Line_T3_C[a].f_t1;
		for (Lbi i=0; i<3; i++) //3t
		{
			Lbd **ftk = &ft[i].F[0];
			for (Lbi j=0; j<Q; j++) //Q
				cudaFree(ftk[j]);
		}
	}

	for (Lbi a=0; a<number_area; a++)
	{
		cudaFree(upsx[a]);
		cudaFree(upsy[a]);

		free(h_uvp2_b[a]);
		cudaFree(d_uvp2_b[a]);

		free(h_uv2_b[a]);
		cudaFree(d_uv2_b[a]);
	}
	/************************************************************************/
	
	cudaFree(d_Point_drag);

	/************************************************************************/

	for (Lbi a=0; a<number_area; a++)
	{
		free(Variables_h[a].solid);
		cudaFree(Variables_d[a].solid);
		Lbd* *Vh = &Variables_h[a].x;
		Lbd* *Vd = &Variables_d[a].x;
		for (Lbi b=0; b<1+3*DMS+2*Q; b++)
		{
			free(Vh[b]);
			Vh[b] = NULL;
			cudaFree(Vd[b]); 
			Vd[b] = NULL;
		}
	}
	free(grid);
	cudaDeviceReset();
	return 0; 
}