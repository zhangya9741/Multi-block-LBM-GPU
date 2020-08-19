#include <omp.h>
#include<stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include <assert.h>
#include <time.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <cusparse.h>
#include <cusparse_v2.h>
#include <string>
using std::string;
//#include "leak_detector_c.h"

#define PI 3.1415926
#define Q 9
#define DMS 2

typedef unsigned int Lbui;
typedef int Lbi;
typedef double Lbd;

#define FLUID 0
#define SOLID_SBB 1
#define SOLID_HBB 3
#define NEUMANN 4
#define BOUNDARY_INLET 5
#define BOUNDARY_OUTLET 7
#define C_BOUNDARY 9
#define F1_BOUNDARY 11
#define F2_BOUNDARY 12
#define IGNORED 13

#define BOUNDARY_NEQ 20
#define ZOUHE_VELOCITY 21
#define ZOUHE_PRESSURE 22
#define ZERO_GRADIENT 23


Lbd w[9]={4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
Lbi e[2*9]={0, 1, 0, -1, 0, 1, -1, -1, 1,0, 0, 1, 0, -1, 1, 1, -1, -1};
Lbi opp[9]={0,3,4,1,2,7,8,5,6};//��׼����
Lbi slip[9]={0,1,4,3,2,8,7,6,5};

__constant__ Lbd w_gpu[Q]={4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
__constant__ Lbi e_gpu[2*Q]={0, 1, 0, -1, 0, 1, -1, -1, 1,0, 0, 1, 0, -1, 1, 1, -1, -1};
__constant__ Lbi opp_gpu[Q]={0,3,4,1,2,7,8,5,6};//��׼����
__constant__ Lbi slip_gpu[Q]={0,1,4,3,2,8,7,6,5};

__constant__ Lbd A_gpu=1.0,B_gpu=-3.0;
__constant__ Lbd M_gpu[Q][Q]={
	{1,1,1,1,1,1,1,1,1},
	{-4,-1,-1,-1,-1,2,2,2,2},
	{4,-2,-2,-2,-2,1,1,1,1},
	{0,1,0,-1,0,1,-1,-1,1},
	{0,-2,0,2,0,1,-1,-1,1},
	{0,0,1,0,-1,1,1,-1,-1},
	{0,0,-2,0,2,1,1,-1,-1},
	{0,1,-1,1,-1,0,0,0,0},
	{0,0,0,0,0,1,-1,1,-1}
};
__constant__ Lbd InvM_gpu[Q][Q]={
	{0.111111111111111,-0.111111111111111,0.111111111111111,0,0,0,0,0,0},
	{0.111111111111111,-0.0277777777777778,-0.0555555555555556,0.166666666666667,-0.166666666666667,0,0,0.250000000000000,0},
	{0.111111111111111,-0.0277777777777778,-0.0555555555555555,0,-1.38777878078145e-17,0.166666666666667,-0.166666666666667,-0.250000000000000,0},
	{0.111111111111111,-0.0277777777777778,-0.0555555555555556,-0.166666666666667,0.166666666666667,0,0,0.250000000000000,0},
	{0.111111111111111,-0.0277777777777778,-0.0555555555555556,0,1.38777878078145e-17,-0.166666666666667,0.166666666666667,-0.250000000000000,0},
	{0.111111111111111,0.0555555555555556,0.0277777777777778,0.166666666666667,0.0833333333333333,0.166666666666667,0.0833333333333333,0,0.250000000000000},
	{0.111111111111111,0.0555555555555556,0.0277777777777778,-0.166666666666667,-0.0833333333333333,0.166666666666667,0.0833333333333333,0,-0.250000000000000},
	{0.111111111111111,0.0555555555555556,0.0277777777777778,-0.166666666666667,-0.0833333333333333,-0.166666666666667,-0.0833333333333333,0,0.250000000000000},
	{0.111111111111111,0.0555555555555556,0.0277777777777778,0.166666666666667,0.0833333333333333,-0.166666666666667,-0.0833333333333333,0,-0.250000000000000},
};