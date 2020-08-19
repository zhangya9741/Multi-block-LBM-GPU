
__global__ void kernel_sum_N_drag(AREA area, D2Q9 d_pp, int *Ns_B)
{
	unsigned int bx=blockIdx.x; unsigned int by=blockIdx.y;
	unsigned int tx=threadIdx.x; unsigned int ty=threadIdx.y;
	int i=bx*blockDim.x+tx;
	int j=by*blockDim.y+ty;

	if (tx == 0 && ty == 0)
		Ns_B[bx+by*gridDim.x] = 0;

	__shared__ int Ns[_BLOCKX*_BLOCKY];
	Ns[tx+ty*_BLOCKX] = 0;

	if (i>=0 && i<area.nx && j>=0 && j<area.ny)
	{
		Lbi ixy = i+j*area.nx;
		Lbi N = 0;
		if (d_pp.solid[ixy] == SOLID_SBB)
		{
			for(Lbi k=0;k<Q;k++)
			{
				Lbi id=i-e_gpu[k];
				Lbi jd=j-e_gpu[k+Q];
				if(id<0 || id>area.nx-1 || jd<0 || jd>area.ny-1)	
					continue;

				Lbi ixyd = id+jd*area.nx;
				if (d_pp.solid[ixyd] == FLUID)
					N++;
			}
		}
		Ns[tx+ty*_BLOCKX] = N;
		__syncthreads();

		for(int ii = _BLOCKX*_BLOCKY / 2; ii > 0; ii >>= 1)
		{
			if(tx+ty*_BLOCKX < ii)
				Ns[tx+ty*_BLOCKX] += Ns[tx+ty*_BLOCKX+ii];
			__syncthreads();
		}

		if(tx == 0 && ty == 0)
		{
			Ns_B[bx+by*gridDim.x] = Ns[0];
		}
	}
}

Lbi sum_N_drag(dim3 block, dim3 grid, AREA area, D2Q9 pp)
{
	int Bx, By;
	Bx = int (ceil((area.nx)/double(_BLOCKX)));
	By = int (ceil((area.ny)/double(_BLOCKY)));

	int *h_Nq = (int *)malloc(Bx*By*sizeof(int));
	int *d_Nq;
	cudaMalloc((void **)&d_Nq, Bx*By*sizeof(int));

	dim3 blocksize(_BLOCKX, _BLOCKY);
	dim3 gridsize(Bx, By, 1);
	kernel_sum_N_drag<<<gridsize, blocksize>>>(area, pp, d_Nq);

	cudaMemcpy(h_Nq, d_Nq, Bx*By*sizeof(int), cudaMemcpyDeviceToHost);

	int N = 0;
	for(int i=0; i<=(Bx*By-1); i++)
	{
		N+=h_Nq[i];
	}
	cudaFree(d_Nq);
	free(h_Nq);

	return N;
}

//fdrag_store保持原状，用cudaMemcpy拷贝
void fdrag_store(AREA area, D2Q9 h_pp, Int3 *d_pdrag, D2Q9 d_pp, Lbi N)
{

	Lbui NN = area.nx*area.ny;
	cudaMemcpy(h_pp.solid, d_pp.solid, NN*sizeof(Lbi), cudaMemcpyDeviceToHost);
	Lbd* *Vh = &h_pp.x;
	Lbd* *Vd = &d_pp.x;
	for (Lbi b=0; b<25; b++)
	{
		cudaMemcpy(Vh[b], Vd[b], NN*sizeof(Lbd), cudaMemcpyDeviceToHost);
	}

	Int3 *pdrag = (Int3 *)malloc(N*sizeof(Int3));
	Lbi ixy, ixyd, ii = 0;
	for (Lbi i=0; i<area.nx; i++)
	{
		for (Lbi j=0; j<area.ny; j++)
		{
			ixy = i+j*area.nx;

			if (h_pp.solid[ixy] == SOLID_SBB)
			{
				for(Lbi k=0;k<Q;k++)
				{
					Lbi id,jd;
					id=i-e[k];
					jd=j-e[k+Q];
					if(id<0|| id>area.nx-1 || jd<0 || jd>area.ny-1)	
						continue;

					ixyd = id+jd*area.nx;
					if (h_pp.solid[ixyd] == FLUID)
					{
						pdrag[ii].i = i;
						pdrag[ii].j = j;
						pdrag[ii].k = k;
						ii++;
					}
				}
			}
			else 
				continue;
		}
	}
	assert(ii == N);

	cudaMemcpy(d_pdrag, pdrag, N*sizeof(Int3), cudaMemcpyHostToDevice);

	free(pdrag);
}

__global__ void kernel_sum_Force(AREA area, D2Q9 d_pp, Int3 *d_pdrag, Lbi N, Lbd *Fx_B, Lbd *Fy_B)
{
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x;
	int ip=bx*blockDim.x+tx;

	__shared__ Lbd Fxs[_BLOCK_1D];
	__shared__ Lbd Fys[_BLOCK_1D];
	if (tx == 0)
	{
		Fx_B[bx] = 0.0;
		Fy_B[bx] = 0.0;
	}
	Fxs[tx] = 0.0;
	Fys[tx] = 0.0;
	if (ip>=0 && ip<N)
	{
		Lbi i = d_pdrag[ip].i;
		Lbi j = d_pdrag[ip].j;
		Lbi k = d_pdrag[ip].k;
		Lbi id = i-e_gpu[k];
		Lbi jd = j-e_gpu[Q+k];
		Lbi kd = opp_gpu[k];
		Lbi ixy = i+j*area.nx;
		Lbi ixyd = id+jd*area.nx;

		Fxs[tx] = (d_pp.f[k][ixy]*e_gpu[k]-d_pp.f[kd][ixyd]*e_gpu[kd]);
		Fys[tx] = (d_pp.f[k][ixy]*e_gpu[k+Q]-d_pp.f[kd][ixyd]*e_gpu[kd+Q]);

		//printf("%d %d %d %f %f\n", i, j, k, d_pp.f[k][ixy], d_pp.f[kd][ixyd]);

		__syncthreads();

		for(int ii = _BLOCK_1D / 2; ii > 0; ii >>= 1)
		{
			if(tx < ii)
			{
				Fxs[tx] += Fxs[tx+ii];
				Fys[tx] += Fys[tx+ii];
			}
			__syncthreads();
		}

		if(tx == 0)
		{
			Fx_B[bx] = Fxs[0];
			Fy_B[bx] = Fys[0];
		}
	}
}

void getDrag_block(Lbi t, Int3 *d_pdrag, AREA area, D2Q9 d_pp, Lbd *force, Lbi N)
{
	Lbi Bx = int (ceil(N/double (_BLOCK_1D)));
	Lbui size_B = Bx*sizeof(Lbd);
	Lbd *h_Fxb = (Lbd *)malloc(size_B);
	Lbd *d_Fxb = NULL;
	cudaMalloc((void **)&d_Fxb, size_B);

	Lbd *h_Fyb = (Lbd *)malloc(size_B);
	Lbd *d_Fyb = NULL;
	cudaMalloc((void **)&d_Fyb, size_B);

	dim3 blocksize(_BLOCK_1D, 1, 1);
	dim3 gridsize(Bx, 1, 1);
	kernel_sum_Force<<<gridsize, blocksize>>>(area, d_pp, d_pdrag, N, d_Fxb, d_Fyb);

	cudaMemcpy(h_Fxb, d_Fxb, size_B, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_Fyb, d_Fyb, size_B, cudaMemcpyDeviceToHost);

	force[0] = 0.0;
	force[1] = 0.0;
	for (Lbi i=0; i<Bx; i++)
	{
		//printf("%d %f %f\n", i, h_Fxb[i], h_Fyb[i]);
		force[0] += h_Fxb[i];
		force[1] += h_Fyb[i];
	}
	cudaFree(d_Fxb);
	cudaFree(d_Fyb);
	free(h_Fxb);
	free(h_Fyb);

	FILE *fp_cdcl;
	fp_cdcl=fopen("fdrag_t.dat","a");
	fprintf(fp_cdcl,"%d\t%g\t%g\t\n", t, force[0]/(rho0*U0_2*U0_2*area.m*R),force[1]/(rho0*U0_2*U0_2*area.m*R));
	fclose(fp_cdcl);
}