
#ifndef _DATA_LBM
#define _DATA_LBM
#include "../lbm/data_type_lbm.h"
#endif

struct FM_spline3{
	Lbi N;
	Lbd *l, *m, *u;
	Lbd *Varibs_M;
	Lbd *rho, *ux, *uy;
	Lbd *F[Q];
};

struct Inf_Ff_store{
	Lbi areac_n;
	Lbi Nc;
	P4 line_C;
	Lbi areaf_n;
	Lbi Nf;
	P4 line_F;
};

struct F_k3{
	Lbd *F[Q];
};
struct Ff_Store{
	Lbi N;
	Lbi marker[3];
	F_k3 f_t1, f_t, f_t_1;
};