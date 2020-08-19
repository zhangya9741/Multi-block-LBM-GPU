
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

const Lbi LEVEL = 2; //总的加密级数
Lbi mi[LEVEL] = {1, 2}; //每级的加密倍数

const Lbi number_area = 2; //总的计算域数量
const Lbi number_line = 1; //粗细网格之间需要信息交换的边界的数量
const Lbi number_line2 = 0;  //同级网格之间需要信息交换的边界的数量*2

Lbi ASP[LEVEL+1] = {0, 1, 2};  //每级的计算域数量
P4d4i AreaZone[number_area] = 
{{0.0,  1.0,  0.0,     0.5,   0, 0, 0, 1},

{0.0,   1.0,  0.5,    1.0,  0, 0, 0, 0}
};  //每块计算域的无量纲起始终止点以及在四个方向上延长的格点数
//（xstart, xend, ystart, yend, left, right, down, up）

Lbi LTCFFCP[LEVEL+1] = {0, 0, 1};//时间插值时Coarse_Fine 信息交换边界数量按级数分布（0, 1）的累加值
Lbi LS3CFFCP[LEVEL+1] = {0, 1, 1}; //空间插值时的分布(1, 0)的累加值
Lbi LCCP[LEVEL+1] = {0, 0, 0}; //同级间(0, 0)的累加值
P3i NeborCF[number_line] = {
	0, 1, 3
}; //三个整数为一组描述需要CF信息交换的线所对应的粗网格域序号，细网格域序号，以及在粗网格上的方位（0-3，对应左右下上）
P3i NeborCC[1] = {}; //三个整数为一组描述需要同级信息交换的线所对应的网格域A序号，网格域B序号，以及在网格上A的方位（0-3，对应左右下上）
//相邻为AB位次交换,此时第三个值实际上没用到
#endif

#ifdef _GEO_2_5

const Lbi LEVEL = 2; //总的加密级数
Lbi mi[LEVEL] = {1, 2}; //每级的加密倍数

const Lbi number_area = 5; //总的计算域数量
const Lbi number_line = 4; //粗细网格之间需要信息交换的边界的数量
const Lbi number_line2 = 8;  //同级网格之间需要信息交换的边界的数量*2

Lbi ASP[LEVEL+1] = {0, 4, 5};  //每级的计算域数量，累加, 4个mi=1的block，1个mi=2的block, --> 0, 4, 1+4

P4d4i AreaZone[number_area] = 
{{0.0,  1.0,  0.75,     1.0,   0, 0, -1, 0}, //The block (mi=1) near the upper boundary, occupying the domain (x/LX, 0-1.0) and (y/LY, 0.75-1.0) 
{0.0,  1.0,  0.0,   0.25,   0, 0, 0, 1}, 
{0.0,   0.25-0.0625,  0.25,    0.75,  0, 1, 0, 0}, 
{0.25+0.0625,   1.0,  0.25,    0.75,  -1, 0, 0, 0},
//The refined block (mi=2) at center, occupying the domain (x/LX, 0.25-0.0625--0.25+0.0625) and (y/LY, 0.25--0.75) 
{0.25-0.0625,   0.25+0.0625,  0.25,    0.75,  0, 0, 0, 0}
};  //每块计算域的无量纲起始终止点以及在四个方向上延长的格点数,低级block向高级block延长，负向延长一个lattice为-1
//（xstart, xend, ystart, yend, left, right, down, up）

Lbi LTCFFCP[LEVEL+1] = {0, 0, 4};//时间插值时Coarse_Fine信息交换边界数量(number_line)按级数分布（0, 4）的累加值--》（0, 0+0, 0+0+4）
Lbi LS3CFFCP[LEVEL+1] = {0, 4, 4}; //空间插值时的分布(4, 0)的累加值--》（0, 0+4, 0+4+0）
Lbi LCCP[LEVEL+1] = {0, 8, 8}; //同级间交换边界数量(number_line2)的累加值（0, 0+8, 0+8+0）
P3i NeborCF[number_line] = {
	0, 4, 2, 
	1, 4, 3, 
	2, 4, 1, 
	3, 4, 0
}; //三个整数为一组描述需要CF信息交换的线所对应的粗网格域序号，细网格域序号，以及在粗网格上的方位（0-3，对应左右下上）
P3i NeborCC[number_line2] = {
	0, 2, 2, 2, 0, 2,
	0, 3, 2, 3, 0, 2,
	1, 2, 3, 2, 1, 3,
	1, 3, 3, 3, 1, 3
}; //三个整数为一组描述需要同级信息交换的线所对应的网格域A序号，网格域B序号，以及在网格上A的方位（0-3，对应左右下上）
//相邻为AB位次交换,此时第三个值实际上没用到
#endif

#ifdef _GEO_4_9

const Lbi LEVEL = 4;    //总的加密级数
Lbi mi[LEVEL] = {1, 2, 2, 2};   //每级的加密倍数

const Lbi number_area = 10;   //总的计算域数量
const Lbi number_line = 17;   //粗细网格之间需要信息交换的边界的数量
const Lbi number_line2 = 12;  //同级网格之间需要信息交换的边界的数量

Lbi ASP[LEVEL+1] = {0, 4, 6, 9, 10};  //每级的计算域数量
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
};   //每块计算域的无量纲起始终止点以及在四个方向上延长的格点数
//（xstart, xend, ystart, yend, left, right, down, up）

Lbi LTCFFCP[LEVEL+1] = {0, 0, 6, 13, 17};  //时间插值时Coarse_Fine 信息交换边界数量按级数分布（0，6，7，4）的累加值
Lbi LS3CFFCP[LEVEL+1] = {0, 10, 14, 17, 17}; //空间插值时的分布(10, 4, 3)的累加值
Lbi LCCP[LEVEL+1] = {0, 8, 8, 12, 12};  //同级间(8, 0, 4, 0)的累加值
P3i NeborCF[number_line] = {
	0, 4, 1, 1, 4, 2, 1, 5, 2, 2, 5, 0, 3, 5, 3, 3, 4, 3,
	1, 6, 2, 1, 7, 2, 3, 7, 3, 3, 8, 3,
	4, 6, 1, 5, 7, 0, 4, 8, 1,
	4, 9, 1,
	6, 9, 2, 7, 9, 0, 8, 9, 3
}; //三个整数为一组描述需要CF信息交换的线所对应的粗网格域序号，细网格域序号，以及在粗网格上的方位（0-3，对应左右下上）
P3i NeborCC[number_line2] = {
	0, 1, 1, 1, 0, 1, 
	2, 1, 0, 1, 2, 0,
	2, 3, 0, 3, 2, 0,
	0, 3, 1, 3, 0, 1,
	7, 6, 0, 6, 7, 0,
	7, 8, 0, 8, 7, 0
}; //三个整数为一组描述需要同级信息交换的线所对应的网格域A序号，网格域B序号，以及在网格上A的方位（0-3，对应左右下上）
//相邻为AB位次交换,此时第三个值实际上没用到
#endif

