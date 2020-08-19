
#ifndef _DATA_MULTI
#define _DATA_MULTI
#include "multi/data_type_multi.h"
#endif

#define _D2Q9

#define _BLOCKX 16
#define _BLOCKY 16
#define _BLOCK_1D (_BLOCKX*_BLOCKY)

#define _BGK

#define _GEO_2_5  
//#define _GEO_4_9  //Fig. 7 of the paper (DOI: 10.1142/S0219876218400029)
//GEO_(Total number of refinement levels)_(Total number of computationsl domains)£¬
//Please define geometry parameters in geo/geometry.h !!!

#define _LOOPS 1 //Total computational time steps = _LOOPS * _T
#define _T 1000 //Output fluid field every _T time steps.

//The Reynolds number
Lbd Re = 10;

//Fluid denstiy
Lbd rho0 = 1.0;

//!!! Please make the parameters *X and *X_gpu have the same value !!!

//The size of computation domain
Lbi LY = 320;
Lbi LX = 768;

//The radius of cylinder
Lbi R = 10;

//The location of cylinder center
Lbi CX = 160;
Lbi CY = 160;

__constant__ Lbi LY_gpu = 320;
__constant__ Lbi LX_gpu = 768;
__constant__ Lbi R_gpu = 10;
__constant__ Lbi CX_gpu = 160;
__constant__ Lbi CY_gpu = 160;

Lbd U0[2] = { 0.1, 0.0 };
Lbd U0_2 = sqrt(U0[0] * U0[0] + U0[1] * U0[1]);

__constant__ Lbd rho0_gpu = 1.0;
__constant__ Lbd U0_gpu[2] = {0.1, 0.0};

/************************************************************************/