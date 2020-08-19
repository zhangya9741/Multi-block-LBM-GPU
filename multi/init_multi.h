

__global__ void init_Fstore(AREA area, D2Q9 pp, Inf_Ff_store line_Inf, Ff_Store line_t3)
{
	Lbui bx = blockIdx.x;
	Lbui tx = threadIdx.x;
	Lbi ip = bx*blockDim.x+tx;
	Lbi nx = area.nx;

	if (ip < line_t3.N)
	{
		Lbi is = line_Inf.line_F.xstart;
		Lbi ie = line_Inf.line_F.xend;
		Lbi js = line_Inf.line_F.ystart;
		Lbi i = 0, j = 0, ixy = 0;

		if (is == ie){ i = is; j = js+ip;}
		else { i = is+ip; j = js;}
		ixy = i+j*nx;

		for (Lbi k=0; k<Q; k++)
		{
			line_t3.f_t1.F[k][ip] = pp.F[k][ixy];
			line_t3.f_t.F[k][ip] = pp.F[k][ixy];
			line_t3.f_t_1.F[k][ip] = pp.F[k][ixy];
		}
	}
}

__global__ void init_tri_A(AREA area, D2Q9 pp, Inf_Ff_store line_Inf, FM_spline3 line_S3)
{
	Lbui bx = blockIdx.x;
	Lbui tx = threadIdx.x;
	Lbi ip = bx*blockDim.x+tx;

	if (ip < line_S3.N)
	{
		line_S3.l[ip] = 0.5;
		line_S3.m[ip] = 2.0;
		line_S3.u[ip] = 0.5;

		if (ip == 0) line_S3.l[ip] = 0.0;
		if (ip == line_Inf.Nc-1) line_S3.u[ip] = 0.0;

		if (ip == 0 || ip == line_Inf.Nc-2)
		{
			line_S3.u[ip] = 0.0;
			line_S3.l[ip+1] = 0.0;
		}
	}
}