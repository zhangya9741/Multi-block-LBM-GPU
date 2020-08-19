

__global__ void error_L2_area(AREA area, Lbd *x, Lbd *y, Lbd *ux, Lbd *uy, Lbd *upx, Lbd *upy, Lbd *uvp2_B, Lbd *uv2_B)
{
	unsigned int bx=blockIdx.x;  unsigned int by=blockIdx.y;
	unsigned int tx=threadIdx.x; unsigned int ty=threadIdx.y;
	Lbui tip = tx+blockDim.x*ty;
	Lbui bip = bx+by*gridDim.x;
	Lbi i = bx*blockDim.x+tx;
	Lbi j = by*blockDim.y+ty;

	__shared__ Lbd uvp2_s[_BLOCKX*_BLOCKY];
	__shared__ Lbd uv2_s[_BLOCKX*_BLOCKY];
	if (tip == 0)
	{
		uvp2_B[bip] = 0.0;
		uv2_B[bip] = 0.0;
	}
	uvp2_s[tip] = 0.0;
	uv2_s[tip] = 0.0;

	if (i>=0 && i<=area.nx-1 && j>=0 && j<=area.ny-1)
	{
		Lbi ixy = i+j*area.nx;
		Lbd ux_tmp, upx_tmp, uy_tmp, upy_tmp;
		ux_tmp = ux[ixy]; 
		upx_tmp = upx[ixy]; 
		uy_tmp = uy[ixy]; 
		upy_tmp = upy[ixy]; 

		uvp2_s[tip] = (ux_tmp-upx_tmp)*(ux_tmp-upx_tmp)+(uy_tmp-upy_tmp)*(uy_tmp-upy_tmp);
		uv2_s[tip] = ux_tmp*ux_tmp+uy_tmp*uy_tmp;

		upx[ixy] = ux_tmp;
		upy[ixy] = uy_tmp;
	}

	__syncthreads();

	for(Lbi ii = _BLOCKX*_BLOCKY / 2; ii > 0; ii >>= 1)
	{
		if(tip < ii)
		{
			uvp2_s[tip] += uvp2_s[tip+ii];
			uv2_s[tip] += uv2_s[tip+ii];
		}
		__syncthreads();
	}

	if(tip == 0)
	{
		uvp2_B[bip] = uvp2_s[0];
		uv2_B[bip] = uv2_s[0];
	}
}

void error_L2(Lbi t, const Lbi N_areas, dim3 block, dim3 *grid, AREA *areas, D2Q9 *pps_d, Lbd **upsx, Lbd **upsy,
	Lbd **h_uvp2_b, Lbd **d_uvp2_b, Lbd **h_uv2_b, Lbd **d_uv2_b)
{
	Lbd uvp2_all, uv2_all, uvp2_area, uv2_area;
	uvp2_all = 0.0;
	uv2_all = 0.0;

	for (Lbi a=0; a<N_areas; a++)
	{
		error_L2_area<<<grid[a], block>>>(areas[a], pps_d[a].x, pps_d[a].y, pps_d[a].ux, pps_d[a].uy, upsx[a], upsy[a], d_uvp2_b[a], d_uv2_b[a]);

		Lbi Bx = int (ceil(areas[a].nx/double (_BLOCKX)));
		Lbi By = int (ceil(areas[a].ny/double (_BLOCKY)));
		Lbui size_B = Bx*By*sizeof(Lbd);
		cudaMemcpy(h_uvp2_b[a], d_uvp2_b[a], size_B, cudaMemcpyDeviceToHost);
		cudaMemcpy(h_uv2_b[a], d_uv2_b[a], size_B, cudaMemcpyDeviceToHost);

		uvp2_area = 0.0;
		uv2_area = 0.0;
		for (Lbi i=0; i<Bx*By; i++)
		{
			uvp2_area += h_uvp2_b[a][i];
			uv2_area += h_uv2_b[a][i];
		}

		uvp2_all += uvp2_area/areas[a].m;
		uv2_all += uv2_area/areas[a].m;
	}
	FILE *fp_error_L2;
	fp_error_L2=fopen("error_L2.dat","a");
	fprintf(fp_error_L2,"%d\t%g\n", t, sqrt(uvp2_all)/sqrt(uv2_all)); //+0.000001));
	fclose(fp_error_L2);

	/************************************************************************/
}