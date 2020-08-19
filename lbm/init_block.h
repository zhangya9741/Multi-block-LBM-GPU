
__global__ void init_General(AREA area, D2Q9 pp, Lbd *upx, Lbd *upy)
{
	Lbui bx = blockIdx.x; Lbui by = blockIdx.y;
	Lbui tx = threadIdx.x; Lbui ty = threadIdx.y;
	Lbi i = bx*blockDim.x+tx;
	Lbi j = by*blockDim.y+ty;

	if (i < area.nx && j < area.ny)
	{
		Lbi ixy = i+j*area.nx;
		Lbd dx = 1.0/area.m;
		pp.x[ixy] = i*dx+area.base_x;
		pp.y[ixy] = j*dx+area.base_y;

		/************************************************************************/
		/* cylinder³õÊ¼»¯       */
		/************************************************************************/
		pp.rho[ixy] = rho0_gpu;
		pp.ux[ixy] = U0_gpu[0]; //4.0*U0_gpu*(LY_gpu-1.5-pp.y[ixy])*(pp.y[ixy]-0.5)/(LY_gpu-2.0)/(LY_gpu-2.0); //
		pp.uy[ixy] = U0_gpu[1];
		upx[ixy] = U0_gpu[0]; //4.0*U0_gpu*(LY_gpu-1.5-pp.y[ixy])*(pp.y[ixy]-0.5)/(LY_gpu-2.0)/(LY_gpu-2.0); //
		upy[ixy] = U0_gpu[1];

		pp.forcex[ixy] = 0.0;
		pp.forcey[ixy] = 0.0;

		if ((pp.x[ixy]-CX_gpu)*(pp.x[ixy]-CX_gpu)+(pp.y[ixy]-CY_gpu)*(pp.y[ixy]-CY_gpu) < R_gpu*R_gpu)
			pp.solid[ixy] = SOLID_SBB;
		else 
			if (pp.y[ixy]==0.0)
			pp.solid[ixy] = FLUID;
		else if (pp.y[ixy]==LY_gpu-1.0)
			pp.solid[ixy] = FLUID;
		else if (pp.x[ixy]==0.0) 
			pp.solid[ixy] = FLUID;
		else if (pp.x[ixy]==LX_gpu-1.0)
			pp.solid[ixy] = FLUID;
		else if (i==0 || j==0 || i==area.nx-1 || j==area.ny-1)
			pp.solid[ixy] = C_BOUNDARY;
		else
			pp.solid[ixy] = FLUID;


		for (Lbi k=0; k<Q; k++)
		{
			pp.F[k][ixy] = feq(k, pp.rho[ixy], pp.ux[ixy], pp.uy[ixy]);
			pp.f[k][ixy] = pp.F[k][ixy];
		}
	}
}