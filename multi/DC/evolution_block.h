

__global__ void evolution_area_stb3(AREA area, D2Q9 pp)
{
	Lbui bx = blockIdx.x; Lbui by = blockIdx.y;
	Lbui tx = threadIdx.x; Lbui ty = threadIdx.y;
	Lbi i = bx*blockDim.x+tx;
	Lbi j = by*blockDim.y+ty;

	if (i>=0 && i<=area.nx-1 && j>=0 && j<=area.ny-1)
	{
		Lbi ixy = i+j*area.nx;
		Lbd dt = 1.0/area.m;

		Lbd ftmp[Q];
		for (Lbi k=0; k<Q; k++)
			ftmp[k] = pp.f[k][ixy];

		if (pp.solid[ixy]==FLUID || pp.solid[ixy]==C_BOUNDARY) 
		{
			Lbd rhoadd, uaddx, uaddy, forcex, forcey;

			rhoadd = pp.rho[ixy];
			uaddx = pp.ux[ixy];
			uaddy = pp.uy[ixy];

			forcex = pp.forcex[ixy];
			forcey = pp.forcey[ixy];

			/************************************************************************/
			/*_BGK*/
			double v_tmp_x, v_tmp_y, Forcek;
			for(Lbi k=0;k<Q;k++)
			{
				v_tmp_x = 3.0*(e_gpu[k] - uaddx)+9.0*(e_gpu[k]*uaddx+e_gpu[k+Q]*uaddy)*e_gpu[k];
				v_tmp_y = 3.0*(e_gpu[k+Q]-uaddy)+9.0*(e_gpu[k]*uaddx+e_gpu[k+Q]*uaddy)*e_gpu[k+Q];
				Forcek = (1.0-0.5*area.rv_tau)*w_gpu[k]*(v_tmp_x*forcex+v_tmp_y*forcey);
				//Forcek = (1.0-dt*0.5*area.rv_tau)*w_gpu[k]*(v_tmp_x*forcex+v_tmp_y*forcey);

				pp.F[k][ixy] = ftmp[k]-(ftmp[k]-feq(k,rhoadd,uaddx,uaddy))*area.rv_tau+dt*Forcek;
			}
		}
		else if (pp.solid[ixy] == SOLID_HBB)
		{
			for(int k=0;k<Q;k++)
				pp.F[k][ixy]=ftmp[k];
		}
	}
}
__global__ void evolution_area_stb1(AREA area, D2Q9 pp)
{
	Lbui bx = blockIdx.x; Lbui by = blockIdx.y;
	Lbui tx = threadIdx.x; Lbui ty = threadIdx.y;
	Lbi i = bx*blockDim.x+tx;
	Lbi j = by*blockDim.y+ty;

	if (i>=0 && i<=area.nx-1 && j>=0 && j<=area.ny-1)
	{
		Lbi ixy = i+j*area.nx;
		Lbd ftmp[Q];
		for (Lbi k=0; k<Q; k++)
			ftmp[k] = 0.0;
		Lbi id, jd, ixyd;

		//stream
		if (pp.solid[ixy]==FLUID || pp.solid[ixy]==C_BOUNDARY) 
		{
			for(Lbi k=0;k<Q;k++)
			{
				id=i-e_gpu[k];
				jd=j-e_gpu[k+Q];
				if(id<0 || id>area.nx-1 || jd<0 || jd>area.ny-1)	
					continue;

				ixyd = id+jd*area.nx;
				if (pp.solid[ixyd] == SOLID_HBB)  //Hbb
				{
					//moving boundary
					Lbd ubdote = pp.ux[ixyd]*e_gpu[opp_gpu[k]]+pp.uy[ixyd]*e_gpu[opp_gpu[k]+Q];
					ftmp[k] = pp.F[opp_gpu[k]][ixy]-6.0*rho0_gpu*w_gpu[opp_gpu[k]]*ubdote;
				}
				else
					ftmp[k] = pp.F[k][ixyd];
			}
		}
		else if (pp.solid[ixy] == SOLID_HBB)
		{
			for(int k=0;k<Q;k++)
				ftmp[k] = pp.F[k][ixy];
		}

		__syncthreads();

		for(int k=0;k<Q;k++)
			pp.f[k][ixy]=ftmp[k];

		/********************************************************************/
	}
}

__global__ void evolution_area_stb2(AREA area, D2Q9 pp)
{
	Lbui bx = blockIdx.x; Lbui by = blockIdx.y;
	Lbui tx = threadIdx.x; Lbui ty = threadIdx.y;
	Lbi i = bx*blockDim.x+tx;
	Lbi j = by*blockDim.y+ty;

	if (i>=0 && i<=area.nx-1 && j>=0 && j<=area.ny-1)
	{
		Lbi ixy = i+j*area.nx;
		Lbd ftmp[Q];
		for (Lbi k=0; k<Q; k++)
			ftmp[k] = pp.f[k][ixy];

		if (pp.solid[ixy] == FLUID) 
		{
			Lbd rhoadd, uaddx, uaddy;

			Macros(rhoadd, uaddx, uaddy, ftmp);

			pp.rho[ixy] = rhoadd;
			pp.ux[ixy] = uaddx;
			pp.uy[ixy] = uaddy;
		}
	}
}

void evolution_area_stb( dim3 block, dim3 grid, AREA area, D2Q9 pp)
{
	evolution_area_stb3<<<grid, block>>>(area, pp);
	evolution_area_stb1<<<grid, block>>>(area, pp);
	evolution_area_stb2<<<grid, block>>>(area, pp);
}