
#ifndef _DATA_MULTI
#define _DATA_MULTI
#include "multi/data_type_multi.h"
#endif

#ifdef _GEO_1_1

const Lbi LEVEL = 1;
Lbi mi[LEVEL] = {1};

const Lbi number_area = 1;
const Lbi number_line = 0;
const Lbi number_line2 = 0;

Lbi ASP[LEVEL+1] = {0, 1};
P4d4i AreaZone[number_area] = 
{{0.0,   1.0,   0,    1.0,   0,  0, 0, 0}
};

Lbi LTCFFCP[LEVEL+1] = {0, 0};
Lbi LS3CFFCP[LEVEL+1] = {0, 0};
Lbi LCCP[LEVEL+1] = {0, 0};
P3i NeborCF[1] = {};
P3i NeborCC[1] = {};

#endif

#ifdef _GEO_2_2

const Lbi LEVEL = 2;
Lbi mi[LEVEL] = {1, 2};

const Lbi number_area = 2; 
const Lbi number_line = 1; 
const Lbi number_line2 = 0;

Lbi ASP[LEVEL+1] = {0, 1, 2}; 
P4d4i AreaZone[number_area] = 
{{0.0,  1.0,  0.0,     0.5,   0, 0, 0, 1},

{0.0,   1.0,  0.5,    1.0,  0, 0, 0, 0}
}; 
//£¨xstart, xend, ystart, yend, left, right, down, up£©

Lbi LTCFFCP[LEVEL+1] = {0, 0, 1};
Lbi LS3CFFCP[LEVEL+1] = {0, 1, 1};
Lbi LCCP[LEVEL+1] = {0, 0, 0}; 
P3i NeborCF[number_line] = {
	0, 1, 3
};
P3i NeborCC[1] = {}; 
#endif

#ifdef _GEO_2_5

const Lbi LEVEL = 2; //Totol refinement levels
Lbi mi[LEVEL] = {1, 2}; //Refinement ratio of each level, level 0, mi=1; level 1, mi=2

const Lbi number_area = 5; //Total number of blocks of MB-LBM
const Lbi number_line = 4; //Number of coarse-fine boundaries between coarse and fine blocks
const Lbi number_line2 = 8;  //2 * Number of boundaries connecting blocks at the same refinement level, such as coarse-coarse and fine-fine

Lbi ASP[LEVEL+1] = {0, 4, 5};  //Number of blocks at different refinement level (4(mi=1), 1(mi=2)) --> (0, 4, 1+4)

P4d4i AreaZone[number_area] = 
{{0.0,  1.0,  0.75,     1.0,   0, 0, -1, 0}, //The block 0 (mi=1) near the upper boundary, occupying the range (x/LX, 0-1.0) and (y/LY, 0.75-1.0) 
{0.0,  1.0,  0.0,   0.25,   0, 0, 0, 1}, //The block 1 (mi=1)
{0.0,   0.25-0.0625,  0.25,    0.75,  0, 1, 0, 0}, //The block 2 (mi=1)
{0.25+0.0625,   1.0,  0.25,    0.75,  -1, 0, 0, 0},//The block 3 (mi=1)
//The refined block 4 (mi=2) at center, occupying the domain (x/LX, 0.25-0.0625--0.25+0.0625) and (y/LY, 0.25--0.75) 
{0.25-0.0625,   0.25+0.0625,  0.25,    0.75,  0, 0, 0, 0}
}; 
//£¨xstart/LX, xend/LX, ystart/LY, yend/LY, left, right, down, up£©
// The relative location of nearby high level block, (left, right, down, up£©
// Extend/shrink (1/-1) the lower level block to its nearby high level block
// For example, extend toward left (xstart/LX, xend/LX, ystart/LY, yend/LY, -1, 0, 0, 0)

Lbi LTCFFCP[LEVEL+1] = {0, 0, 4};
//Number of boundaries that need temporal interpolation at different level (0, 4) --> (0, 0+0, 0+0+4), 
//total number of these boundaries is (number_line)
Lbi LS3CFFCP[LEVEL+1] = {0, 4, 4}; 
//Number of boundaries that need spatial interpolation at different level (4, 0) --> (0, 0+4, 0+4+0), 
//total number of these boundaries is (number_line)
Lbi LCCP[LEVEL+1] = {0, 8, 8}; 
//Number of boundaries connecting the blocks at the same level (8, 0) --> (0, 8+0, 0+8+0)
//total number of these boundaries is (number_line2)
P3i NeborCF[number_line] = {
	0, 4, 2, //Coarse block ID, Fine block ID, the relative location of the fine block[*(left, right, down, up£©--> (0, 1, 2, 3)]
	1, 4, 3, 
	2, 4, 1, 
	3, 4, 0
}; //Coarse-fine boundaries and their neighbouring blocks
P3i NeborCC[number_line2] = { 
	0, 2, 2, 2, 0, 2, //Block A, Block B, the relative location of Block B[*(left, right, down, up£©--> (0, 1, 2, 3)], Block B, Block A, the relative location of Block A [not used in practice]
	0, 3, 2, 3, 0, 2,
	1, 2, 3, 2, 1, 3,
	1, 3, 3, 3, 1, 3
}; //Boundaries connecting blocks at the same level
#endif

#ifdef _GEO_4_9

const Lbi LEVEL = 4;   
Lbi mi[LEVEL] = {1, 2, 2, 2};  

const Lbi number_area = 10;  
const Lbi number_line = 17; 
const Lbi number_line2 = 12; 

Lbi ASP[LEVEL+1] = {0, 4, 6, 9, 10};  
P4d4i AreaZone[number_area] = 
{{0.0,  0.1,  0,     1.0,   0,  1, 0, 0},
{0.1,   0.4,  0.75,  1.0,   0, 0, -1, 0}, 
{0.4,   1.0,  0.0,   1.0,   -1, 0, 0, 0}, 
{0.1,   0.4,  0.0,   0.25,  0, 0, 0,  1}, 

{0.1,   0.15, 0.25,  0.75,  0,  1, 0, 0},
{0.35,  0.4,  0.25,  0.75,  -1, 0, 0, 0},

{0.15,  0.25, 0.625, 0.75,  0, 0, -1, 0},
{0.25,  0.35, 0.25,  0.75,  -1, 0, 0, 0},
{0.15,  0.25, 0.25,  0.375, 0, 0, 0, 1},

{0.15,  0.25, 0.375, 0.625, 0,  0, 0, 0}
};   
//£¨xstart, xend, ystart, yend, left, right, down, up£©

Lbi LTCFFCP[LEVEL+1] = {0, 0, 6, 13, 17};  
Lbi LS3CFFCP[LEVEL+1] = {0, 10, 14, 17, 17};
Lbi LCCP[LEVEL+1] = {0, 8, 8, 12, 12}; 
P3i NeborCF[number_line] = {
	0, 4, 1, 1, 4, 2, 1, 5, 2, 2, 5, 0, 3, 5, 3, 3, 4, 3,
	1, 6, 2, 1, 7, 2, 3, 7, 3, 3, 8, 3,
	4, 6, 1, 5, 7, 0, 4, 8, 1,
	4, 9, 1,
	6, 9, 2, 7, 9, 0, 8, 9, 3
}; 
P3i NeborCC[number_line2] = {
	0, 1, 1, 1, 0, 1, 
	2, 1, 0, 1, 2, 0,
	2, 3, 0, 3, 2, 0,
	0, 3, 1, 3, 0, 1,
	7, 6, 0, 6, 7, 0,
	7, 8, 0, 8, 7, 0
}; 
#endif

