
#ifndef _LBM
#define _LBM
#include "lbm.h"
#endif

struct D2Q9{
	Lbi *solid;
	Lbd *x, *y;
	Lbd *rho, *ux, *uy;
	Lbd *forcex, *forcey;
	Lbd *f[Q];
	Lbd *F[Q];
};
struct P4{
	Lbi xstart;
	Lbi xend;
	Lbi ystart;
	Lbi yend;
};
struct P3i{
	Lbi nc;
	Lbi nf;
	Lbi direction;
};
struct P4d{
	Lbd xstart;
	Lbd xend;
	Lbd ystart;
	Lbd yend;
};
struct P4d4i
{
	P4d reduced;  
	P4 move;
};
struct Int3{
	Lbi i;
	Lbi j;
	Lbi k;
};

struct AREA{
	Lbi m;
	Lbi nx;
	Lbi ny;
	Lbd base_x;
	Lbd base_y;
	Lbd tau;
	Lbd rv_tau;
};
