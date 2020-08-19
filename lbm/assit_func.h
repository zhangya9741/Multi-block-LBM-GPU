
__device__ inline Lbd feq(Lbi k, Lbd rh, Lbd uxx, Lbd uyy)
{
	Lbd eu,uv,feq0;
	eu=e_gpu[k]*uxx+e_gpu[k+Q]*uyy;
	uv=uxx*uxx+uyy*uyy;
	feq0=w_gpu[k]*rh*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	return feq0;
}
__device__ inline void Meq(Lbd &rh, Lbd &uxx, Lbd &uyy, Lbd *meq)
{
	double uv = uxx*uxx+uyy*uyy;
	meq[0]=rh;
	meq[1]=rh*(-2.0+3.0*uv);
	meq[2]=rh*(A_gpu+B_gpu*uv);
	meq[3]=rh*(uxx);
	meq[4]=-meq[3];
	meq[5]=rh*(uyy);
	meq[6]=-meq[5];
	meq[7]=rh*(uxx*uxx-uyy*uyy);
	meq[8]=rh*(uxx*uyy);
}
__device__ inline void Macros(Lbd &rh, Lbd &uxx, Lbd &uyy, Lbd *ftmp)
{
	rh = ftmp[0]+ftmp[1]+ftmp[2]+ftmp[3]+ftmp[4]+ftmp[5]+ftmp[6]+ftmp[7]+ftmp[8];
	uxx = ftmp[1]-ftmp[3]+ftmp[5]-ftmp[6]-ftmp[7]+ftmp[8];
	uyy = ftmp[2]-ftmp[4]+ftmp[5]+ftmp[6]-ftmp[7]-ftmp[8];

	uxx = uxx/rh;
	uyy = uyy/rh;
}