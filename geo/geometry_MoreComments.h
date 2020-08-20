
#ifndef _DATA_MULTI
#define _DATA_MULTI
#include "multi/data_type_multi.h"
#endif


#ifdef _GEO_2_5

const Lbi LEVEL = 2; // total level of refinement//总的加密级数
Lbi mi[LEVEL] = {1, 2};// ratio of fine grids to coarse grids for each level

const Lbi number_area = 5;  // total number of computing area 
const Lbi number_line = 4; // number of boundary for exchangeing data 
                           // from both side of coarse and fine grids
const Lbi number_line2 = 8;// (need *2) number of boundarys for exchanging 
                           // data from both sides of the same refinement level

Lbi ASP[LEVEL+1] = {0, 4, 5};  // computing areas of each refinement level 
                               //explanation for this:
                               // 0 - refinement level  0
                               // 4 - refinement level  0 + 4 = 4
                               // 1 - refinement level  0 + 4 + 1 = 5

P4d4i AreaZone[number_area] = 
{{0.0,  1.0,  0.75,     1.0,   0, 0, -1, 0},// The block (mi=1) near the upper 
                                            // boundary, occupying the domain 
                                            // (x/LX, 0-1.0) & (y/LY, 0.75-1.0) 
{0.0,  1.0,  0.0,   0.25,   0, 0, 0, 1}, 
{0.0,   0.25-0.0625,  0.25,    0.75,  0, 1, 0, 0}, 
{0.25+0.0625,   1.0,  0.25,    0.75,  -1, 0, 0, 0},
// The refined block (mi=2) at center, occupying the domain 
// (x/LX, 0.25-0.0625--0.25+0.0625)  &  (y/LY, 0.25--0.75) 
//
{0.25-0.0625,   0.25+0.0625,  0.25,    0.75,  0, 0, 0, 0}
}; 
// dimensionless start and end points of each computing area,
// and number of points in each four direction extended line. 
//（xstart, xend, ystart, yend, left, right, down, up）
// ************************************************************
// the low-level block is extended to the high-level in both dircetion,
// negative direction = -1  &  positive direction = lattice number + 1


Lbi LTCFFCP[LEVEL+1] = {0, 0, 1};// number of boundarys for exchanging 
                                 // coarse-fine grid sata when 
                                 // do temporal interpolation, 
                                 // （0, 1）

Lbi LTCFFCP[LEVEL+1] = {0, 0, 4};// number of boundarys for exchanging 
                                 // coarse-fine grid sata when 
                                 // do temporal interpolation, 
                                 // cumulative value of each level 
                                 // (0, 0 + 0, 0 + 0 + 4)
                                   
Lbi LS3CFFCP[LEVEL+1] = {0, 4, 4};// number of boundarys for exchanging  
                                  // coarse-fine grid sata when 
                                  // do space interpolation 
                                  // cumulative value of each level 
                                  // (0, 0 + 0, 0 + 0 + 4)
                                    
Lbi LCCP[LEVEL+1] = {0, 8, 8}; // cimilative value of the same level
                               // for exchanging data   (number_line2)
                               // cumulative value of each level 
                               // (0, 0 + 8, 0 + 8 + 0)

P3i NeborCF[number_line] = {
	0, 4, 2, 
	1, 4, 3, 
	2, 4, 1, 
	3, 4, 0
}; // a integer array to discribe grids order which need to be refined. 
   // each array contains three integers, eg:
   // 1) order of coarse grids, 
   // 2) order of fine grids, 
   // 3) direction of coares grids to fine grids (0~3 for LRDU)

P3i NeborCC[number_line2] = {
	0, 2, 2, 2, 0, 2,
	0, 3, 2, 3, 0, 2,
	1, 2, 3, 2, 1, 3,
	1, 3, 3, 3, 1, 3
}; // a integer array to discribe grids order which need to exchanging data.
   // this is exchanging data from the same refinement level, 
   // !!! do not need interpolation 
   // each array contains three integers, eg:
   // 1) order of A grids area, 
   // 2) order of B grids area, 
   // 3) direction of B area to A area (0~3 for LeftRigntDownUp)
   // if two areas is neighbors, just exhanche order of A and B, 
   // do not need the third integer
#endif

