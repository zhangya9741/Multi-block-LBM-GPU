
__global__ void c_c_boundary_transfer(AREA area1, AREA area2, D2Q9 pp1, D2Q9 pp2, Inf_Ff_store line_Inf)
{
	Lbui bx = blockIdx.x;
	Lbui tx = threadIdx.x;
	Lbi ip = bx*blockDim.x+tx;

	Lbi xs1 = line_Inf.line_C.xstart;
	Lbi xe1 = line_Inf.line_C.xend;
	Lbi ys1 = line_Inf.line_C.ystart;
	Lbi ye1 = line_Inf.line_C.yend;
	Lbi nx1 = area1.nx;

	Lbi xs2 = line_Inf.line_F.xstart;
	Lbi ys2 = line_Inf.line_F.ystart;
	Lbi nx2 = area2.nx;

	//1 c为源数据，2 f为目标
	if (ip < line_Inf.Nc)
	{
		Lbi i1 = 0, j1 = 0,  ixy1 = 0, i2 = 0, j2 = 0, ixy2 = 0;
		if (xs1 == xe1)
		{
			i1 = xs1;
			j1 = ys1+ip;
			i2 = xs2;
			j2 = ys2+ip;
		}
		if (ys1 == ye1)
		{
			i1 = xs1+ip;
			j1 = ys1;
			i2 = xs2+ip;
			j2 = ys2;
		}
		ixy1 = i1+j1*nx1;
		ixy2 = i2+j2*nx2;

		pp2.rho[ixy2] = pp1.rho[ixy1];
		pp2.ux[ixy2] = pp1.ux[ixy1];
		pp2.uy[ixy2] = pp1.uy[ixy1];

		for (Lbi k=0; k<Q; k++)
			pp2.f[k][ixy2] = pp1.f[k][ixy1];
	}
}

__global__ void init_Spline3(AREA area, D2Q9 pp, Inf_Ff_store line_Inf, FM_spline3 line_S3)
{	
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x; 
	int ip=bx*blockDim.x+tx;

	if (ip >= 0 && ip < line_S3.N)
	{
		Lbi i, j, ixy;
		Lbi xs = line_Inf.line_C.xstart;
		Lbi xe = line_Inf.line_C.xend;
		Lbi ys = line_Inf.line_C.ystart;
		Lbi nx = area.nx;

		if (xs == xe) 
		{ 
			i = xs; j = ys+ip;
			ixy = i+j*nx;
		}
		else 
		{ 
			i = xs+ip; j = ys;
			ixy = i+j*nx;
		}

		line_S3.rho[ip] = pp.rho[ixy];
		line_S3.ux[ip] = pp.ux[ixy];
		line_S3.uy[ip] = pp.uy[ixy];

		for (Lbi k=0; k<Q; k++)
			line_S3.F[k][ip] = pp.f[k][ixy];

		__syncthreads();

		if (ip == 0 || ip == line_S3.N-1)
		{
			line_S3.Varibs_M[ip] = 0.0;
			line_S3.Varibs_M[line_S3.N+ip] = 0.0;
			line_S3.Varibs_M[2*line_S3.N+ip] = 0.0;

			for (Lbi k=0; k<Q; k++)
				line_S3.Varibs_M[(3+k)*line_S3.N+ip] = 0.0;
		}
		else
		{
			//line_S3.Varibs_M[ip] = 1.5*(line_S3.rho[ip+1]-line_S3.rho[ip-1]);
			//line_S3.Varibs_M[line_S3.N+ip] = 1.5*(line_S3.ux[ip+1]-line_S3.ux[ip-1]);
			//line_S3.Varibs_M[2*line_S3.N+ip] = 1.5*(line_S3.uy[ip+1]-line_S3.uy[ip-1]);

			//for (Lbi k=0; k<Q; k++)
			//	line_S3.Varibs_M[(3+k)*line_S3.N+ip] = 1.5*(line_S3.F[k][ip+1]-line_S3.F[k][ip-1]);

			line_S3.Varibs_M[ip] = (6.0*line_S3.rho[ip]-3.0*line_S3.rho[ip-1]-3.0*line_S3.rho[ip+1]);
			line_S3.Varibs_M[line_S3.N+ip] = (6.0*line_S3.ux[ip]-3.0*line_S3.ux[ip-1]-3.0*line_S3.ux[ip+1]);
			line_S3.Varibs_M[2*line_S3.N+ip] = (6.0*line_S3.uy[ip]-3.0*line_S3.uy[ip-1]-3.0*line_S3.uy[ip+1]);

			for (Lbi k=0; k<Q; k++)
			{
				line_S3.Varibs_M[(3+k)*line_S3.N+ip] = (6.0*line_S3.F[k][ip]-3.0*line_S3.F[k][ip-1]-3.0*line_S3.F[k][ip+1]);
			}
		}
		//printf("%d %f %f %f %f %f %f %f %f %f\n", ip, pp.f[0][ixy], pp.f[1][ixy], pp.f[2][ixy], pp.f[3][ixy], pp.f[4][ixy], pp.f[5][ixy]
		//, pp.f[6][ixy], pp.f[7][ixy], pp.f[8][ixy]);

		//printf("%d %f %f %f %f %f %f %f %f %f\n", ip, line_S3.Varibs_M[(3+0)*line_S3.N+ip], line_S3.Varibs_M[(3+1)*line_S3.N+ip]
		//, line_S3.Varibs_M[(3+2)*line_S3.N+ip], line_S3.Varibs_M[(3+3)*line_S3.N+ip]
		//, line_S3.Varibs_M[(3+4)*line_S3.N+ip], line_S3.Varibs_M[(3+5)*line_S3.N+ip]
		//, line_S3.Varibs_M[(3+6)*line_S3.N+ip], line_S3.Varibs_M[(3+7)*line_S3.N+ip]
		//, line_S3.Varibs_M[(3+8)*line_S3.N+ip]);
	}
}

__device__ inline double Spline3_S(double xx, double x_1, double h, double fi_1, double fi, double mi_1, double mi)
{
	double x, a, b, c, d;
	x = x_1+h;

	//a = mi_1*(xx-x)*(xx-x)*(xx-x_1);
	//b = mi*(xx-x_1)*(xx-x_1)*(xx-x);
	//c = fi_1*(xx-x)*(xx-x)*(1.0+2*(xx-x_1));
	//d = fi*(xx-x_1)*(xx-x_1)*(1.0+2*(x-xx));

	a = mi_1*(x-xx)*(x-xx)*(x-xx)*0.1666667;
	b = mi*(xx-x_1)*(xx-x_1)*(xx-x_1)*0.1666667;
	c = (fi_1-mi_1*h*h*0.1666667)*(x-xx);
	d = (fi-mi*h*h*0.1666667)*(xx-x_1);

	return (a+b+c+d);
}

__global__ void calculate_Spline3_gpu(AREA areac, AREA areaf, D2Q9 pp, 
	Inf_Ff_store line_Inf, F_k3 lt3, FM_spline3 line_S3)
{
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x; 
	int ip=bx*blockDim.x+tx;

	if (ip < line_Inf.Nf)
	{
		Lbi i, j, ixy;
		Lbi xs = line_Inf.line_F.xstart;
		Lbi xe = line_Inf.line_F.xend;
		Lbi ys = line_Inf.line_F.ystart;

		if (xs == xe) 
		{ 
			i = xs; j = ys+ip;
			ixy = i+j*areaf.nx;
		}
		else 
		{ 
			i = xs+ip; j = ys;
			ixy = i+j*areaf.nx;
		}
		if (pp.solid[ixy] != SOLID_HBB)
		{
			Lbd xx = Lbd (ip)/Lbd (areaf.m/areac.m);
			Lbi I = ceil(xx);

			if (ceil(xx) == 0)
				I = Lbi (ceil(xx)+1.0);

			Lbd rho_tmp = Spline3_S(xx, (I-1.0), 1.0, line_S3.rho[I-1],  line_S3.rho[I], line_S3.Varibs_M[I-1],  line_S3.Varibs_M[I]);
			Lbd ux_tmp = Spline3_S(xx, (I-1.0), 1.0, line_S3.ux[I-1],  line_S3.ux[I], line_S3.Varibs_M[line_S3.N+I-1],  line_S3.Varibs_M[line_S3.N+I]);
			Lbd uy_tmp = Spline3_S(xx, (I-1.0), 1.0, line_S3.uy[I-1],  line_S3.uy[I], line_S3.Varibs_M[2*line_S3.N+I-1],  line_S3.Varibs_M[2*line_S3.N+I]);

			//printf("%d %f %f %f\n", ip, rho_tmp, ux_tmp, uy_tmp);

			//Lbd para = (areaf.tau-1.0)/(areaf.m/areac.m*(areac.tau-1.0));
			Lbd para = areaf.tau/(areaf.m/areac.m*areac.tau);
			Lbd f_tmp, feq_tmp;

			for (Lbi k=0; k<Q; k++)
			{
				f_tmp = Spline3_S(xx, (I-1.0), 1.0, line_S3.F[k][I-1], line_S3.F[k][I], line_S3.Varibs_M[(3+k)*line_S3.N+I-1], line_S3.Varibs_M[(3+k)*line_S3.N+I]);
				feq_tmp = feq(k, rho_tmp, ux_tmp, uy_tmp);
				lt3.F[k][ip] = feq_tmp+para*(f_tmp-feq_tmp);
			}
		}
	}
}

void space_interpolation_Spline3(AREA areac, AREA areaf, D2Q9 ppc, D2Q9 ppf,
	Inf_Ff_store line_Inf, Ff_Store line_t3, FM_spline3 line_S3, cusparseHandle_t handle)
{
	init_Spline3<<<(int (ceil(line_Inf.Nc/double (_BLOCK_1D)))), _BLOCK_1D>>>
		(areac, ppc, line_Inf, line_S3);

	cusparseDgtsv(handle, line_Inf.Nc, 12, line_S3.l, line_S3.m, line_S3.u, line_S3.Varibs_M, line_Inf.Nc);

	F_k3 *lt3t = &line_t3.f_t1;
	F_k3 lt3t_t = lt3t[line_t3.marker[0]];
	calculate_Spline3_gpu<<<(int (ceil(line_Inf.Nf/double (_BLOCK_1D)))), _BLOCK_1D>>>
		(areac, areaf, ppf, line_Inf, lt3t_t, line_S3);
}

__global__ void temporal_interpolation(AREA area, D2Q9 pp, Inf_Ff_store line_Inf, Lbi N,
	F_k3 Ft_1, F_k3 Ft, F_k3 Ft1, Lbd xm)
{
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x; 
	int ip=bx*blockDim.x+tx;
	//if (ip == 0) printf("%d %f\n", ip, xm);

	if (ip < N)
	{
		Lbi i, j;
		Lbi xs = line_Inf.line_F.xstart;
		Lbi xe = line_Inf.line_F.xend;
		Lbi ys = line_Inf.line_F.ystart;

		if (xs == xe) { i = xs; j = ys+ip;}
		else { i = xs+ip; j = ys;}
		Lbi ixy = i+j*area.nx;
		if (pp.solid[ixy] != SOLID_HBB)
		{
			Lbd a, b, c;
			a = 0.5*(xm-1.0)*xm;
			b = (xm-1.0)*(xm+1.0);
			c = 0.5*xm*(xm+1.0);

			pp.f[0][ixy] = a*Ft_1.F[0][ip]-b*Ft.F[0][ip]+c*Ft1.F[0][ip];
			pp.f[1][ixy] = a*Ft_1.F[1][ip]-b*Ft.F[1][ip]+c*Ft1.F[1][ip];
			pp.f[2][ixy] = a*Ft_1.F[2][ip]-b*Ft.F[2][ip]+c*Ft1.F[2][ip];
			pp.f[3][ixy] = a*Ft_1.F[3][ip]-b*Ft.F[3][ip]+c*Ft1.F[3][ip];
			pp.f[4][ixy] = a*Ft_1.F[4][ip]-b*Ft.F[4][ip]+c*Ft1.F[4][ip];
			pp.f[5][ixy] = a*Ft_1.F[5][ip]-b*Ft.F[5][ip]+c*Ft1.F[5][ip];
			pp.f[6][ixy] = a*Ft_1.F[6][ip]-b*Ft.F[6][ip]+c*Ft1.F[6][ip];
			pp.f[7][ixy] = a*Ft_1.F[7][ip]-b*Ft.F[7][ip]+c*Ft1.F[7][ip];
			pp.f[8][ixy] = a*Ft_1.F[8][ip]-b*Ft.F[8][ip]+c*Ft1.F[8][ip];

			Lbd rho_tmp=0.0, ux_tmp=0.0, uy_tmp=0.0;
			for (Lbi k=0; k<Q; k++)
			{
				rho_tmp += pp.f[k][ixy]; 
				ux_tmp += pp.f[k][ixy]*e_gpu[k]; 
				uy_tmp += pp.f[k][ixy]*e_gpu[k+Q]; 
			}
			pp.rho[ixy] = rho_tmp;
			pp.ux[ixy] = ux_tmp;
			pp.uy[ixy] = uy_tmp;
		}
	}
}

__global__ void c_f_boundary_transfer(AREA area, D2Q9 pp, Inf_Ff_store line_Inf, Lbi N, F_k3 Ft)
{
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x; 
	int ip=bx*blockDim.x+tx;

	if (ip < N)
	{
		Lbi i, j;
		Lbi xs = line_Inf.line_F.xstart;
		Lbi xe = line_Inf.line_F.xend;
		Lbi ys = line_Inf.line_F.ystart;

		if (xs == xe) { i = xs; j = ys+ip;}
		else { i = xs+ip; j = ys;}
		Lbi ixy = i+j*area.nx;
		if (pp.solid[ixy] != SOLID_HBB)
		{
			for (Lbi k=0; k<Q; k++)
			{
				pp.f[k][ixy] = Ft.F[k][ip];
			}

			Lbd rho_tmp=0.0, ux_tmp=0.0, uy_tmp=0.0;
			for (Lbi k=0; k<Q; k++)
			{
				rho_tmp += Ft.F[k][ip]; 
				ux_tmp += Ft.F[k][ip]*e_gpu[k]; 
				uy_tmp += Ft.F[k][ip]*e_gpu[k+Q]; 
			}
			pp.rho[ixy] = rho_tmp;
			pp.ux[ixy] = ux_tmp;
			pp.uy[ixy] = uy_tmp;
		}
	}
}

__global__ void f_c_boundary_transfer(AREA areac, AREA areaf, D2Q9 ppc, D2Q9 ppf, Inf_Ff_store line_Inf)
{
	unsigned int bx=blockIdx.x;
	unsigned int tx=threadIdx.x; 
	int ip=bx*blockDim.x+tx;

	if (ip < line_Inf.Nc)
	{
		Lbi i_c, j_c, i_f, j_f;
		Lbi xs_c = line_Inf.line_C.xstart;
		Lbi xe_c = line_Inf.line_C.xend;
		Lbi ys_c = line_Inf.line_C.ystart;

		Lbi xs_f = line_Inf.line_F.xstart;
		Lbi ys_f = line_Inf.line_F.ystart;

		if (xs_c == xe_c) {
			i_c = xs_c; j_c = ys_c+ip;
			i_f = xs_f; j_f = ys_f+areaf.m/areac.m*ip;
		}
		else { 
			i_c = xs_c+ip; j_c = ys_c;
			i_f = xs_f+areaf.m/areac.m*ip; j_f = ys_f;
		}
		Lbi ixyc = i_c+j_c*areac.nx;
		Lbi ixyf = i_f+j_f*areaf.nx;
		if (ppc.solid[ixyc] != SOLID_HBB)
		{
			ppc.rho[ixyc] = ppf.rho[ixyf];
			ppc.ux[ixyc] = ppf.ux[ixyf];
			ppc.uy[ixyc] = ppf.uy[ixyf];

			//Lbd para = areaf.m/areac.m*(areac.tau-1.0)/(areaf.tau-1.0);
			Lbd para = areaf.m/areac.m*areac.tau/areaf.tau;
			Lbd f_tmp, feq_tmp;
			for (Lbi k=0; k<Q; k++)
			{
				f_tmp = ppf.F[k][ixyf];
				feq_tmp = feq(k, ppf.rho[ixyf], ppf.ux[ixyf], ppf.uy[ixyf]);
				ppc.f[k][ixyc] = feq_tmp+para*(f_tmp-feq_tmp);
			}
		}
	}
}