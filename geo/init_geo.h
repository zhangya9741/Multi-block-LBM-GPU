

#ifndef _GEO
#define _GEO
#include "geometry.h"
#endif

// Geometry arrangement    

P4 AreaDefine[number_area];
AREA Areas[number_area];
Inf_Ff_store Line_Inf_C[(number_line<=0) ? 1 : number_line];
Inf_Ff_store Line_Inf_F[(number_line<=0) ? 1 : number_line];
Inf_Ff_store Line_Inf_CC[(number_line2<=0) ? 1 : number_line2]; 


void Line_init(Lbi Direction, P4 &areaD_C, AREA &areas_C, Inf_Ff_store &lines_C, P4 &areaD_F, AREA &areas_F, Inf_Ff_store &lines_F)
{

	if (Direction == 0 || Direction == 1)
	{
		lines_C.line_C.xstart = (Direction==0 ? areas_C.m : areas_C.nx-1-areas_C.m);
		lines_C.line_C.xend = lines_C.line_C.xstart;

		lines_C.line_F.xstart = (Direction==0 ? areas_F.nx-1 : 0);
		lines_C.line_F.xend = lines_C.line_F.xstart;

		lines_F.line_C.xstart = (Direction==0 ? 0 : areas_C.nx-1); 
		lines_F.line_C.xend = lines_F.line_C.xstart;

		lines_F.line_F.xstart = (Direction==0 ? areas_F.nx-areas_F.m-1 : areas_F.m);
		lines_F.line_F.xend = lines_F.line_F.xstart;

		lines_C.line_C.ystart = (areaD_F.ystart-areaD_C.ystart<=0) ? 0 : (int)(areaD_F.ystart-areaD_C.ystart)*areas_C.m;
		lines_C.line_C.yend = (areaD_F.yend-areaD_C.yend>=0) ? areas_C.ny-1 : (int)(areaD_F.yend-areaD_C.ystart)*areas_C.m;

		lines_C.line_F.ystart = (areaD_C.ystart-areaD_F.ystart<=0) ? 0 : (int)(areaD_C.ystart-areaD_F.ystart)*areas_F.m;
		lines_C.line_F.yend = (areaD_C.yend-areaD_F.yend>=0) ? areas_F.ny-1 : (int)(areaD_C.yend-areaD_F.ystart)*areas_F.m;

		lines_F.line_C.ystart = lines_C.line_C.ystart;
		lines_F.line_C.yend = lines_C.line_C.yend;

		lines_F.line_F.ystart = lines_C.line_F.ystart;
		lines_F.line_F.yend = lines_C.line_F.yend;

		lines_C.Nc = lines_C.line_C.yend-lines_C.line_C.ystart+1;
		lines_C.Nf = lines_C.line_F.yend-lines_C.line_F.ystart+1;

		lines_F.Nc = lines_C.Nc;
		lines_F.Nf = lines_C.Nf;
	}

	else if (Direction == 2 || Direction == 3)
	{

		lines_C.line_C.ystart = (Direction==2 ? areas_C.m : areas_C.ny-1-areas_C.m); 
		lines_C.line_C.yend = lines_C.line_C.ystart;

		lines_C.line_F.ystart = (Direction==2 ? areas_F.ny-1 : 0); 
		lines_C.line_F.yend = lines_C.line_F.ystart;

		lines_F.line_C.ystart = (Direction==2 ? 0 : areas_C.ny-1); 
		lines_F.line_C.yend = lines_F.line_C.ystart;

		lines_F.line_F.ystart = (Direction==2 ? areas_F.ny-areas_F.m-1 : areas_F.m); 
		lines_F.line_F.yend = lines_F.line_F.ystart;

		lines_C.line_C.xstart = (areaD_F.xstart-areaD_C.xstart<=0) ? 0 : (int)(areaD_F.xstart-areaD_C.xstart)*areas_C.m;
		lines_C.line_C.xend = (areaD_F.xend-areaD_C.xend>=0) ? areas_C.nx-1 : (int)(areaD_F.xend-areaD_C.xstart)*areas_C.m;

		lines_C.line_F.xstart = (areaD_C.xstart-areaD_F.xstart<=0) ? 0 : (int)(areaD_C.xstart-areaD_F.xstart)*areas_F.m;
		lines_C.line_F.xend = (areaD_C.xend-areaD_F.xend>=0) ? areas_F.nx-1 : (int)(areaD_C.xend-areaD_F.xstart)*areas_F.m;

		lines_F.line_C.xstart = lines_C.line_C.xstart;
		lines_F.line_C.xend = lines_C.line_C.xend;

		lines_F.line_F.xstart = lines_C.line_F.xstart;
		lines_F.line_F.xend = lines_C.line_F.xend;

		lines_C.Nc = lines_C.line_C.xend-lines_C.line_C.xstart+1;
		lines_C.Nf = lines_C.line_F.xend-lines_C.line_F.xstart+1;

		lines_F.Nc = lines_C.Nc;
		lines_F.Nf = lines_C.Nf;
	}
}


void geo_init()
{
	for (Lbi i=0; i<number_area; i++)
	{
		AreaDefine[i].xstart = (int)(AreaZone[i].reduced.xstart*(LX-1))+AreaZone[i].move.xstart;
		AreaDefine[i].xend = (int)(AreaZone[i].reduced.xend*(LX-1))+AreaZone[i].move.xend;
		AreaDefine[i].ystart = (int)(AreaZone[i].reduced.ystart*(LY-1))+AreaZone[i].move.ystart;
		AreaDefine[i].yend = (int)(AreaZone[i].reduced.yend*(LY-1))+AreaZone[i].move.yend;
	}

	Lbi mm = 1;

	Lbd niu = U0_2*2*R/Re;//U0_2*(LX-2)/Re;

	Lbd tau = 0.5+3.0*niu;
	Lbd rv_tau = 1.0/tau;
	for (Lbi i=0; i<LEVEL; i++)
	{
		mm *= mi[i];
		tau = 0.5+mi[i]*(tau-0.5);
		rv_tau = 1.0/tau;
		for (Lbi j=ASP[i]; j<ASP[i+1]; j++)
		{
			Areas[j].m = mm;
			Areas[j].nx = (Lbi)(AreaDefine[j].xend-AreaDefine[j].xstart)*mm+1;
			Areas[j].ny = (Lbi)(AreaDefine[j].yend-AreaDefine[j].ystart)*mm+1;
			Areas[j].base_x = AreaDefine[j].xstart;
			Areas[j].base_y = AreaDefine[j].ystart;
			Areas[j].tau = tau;
			Areas[j].rv_tau = rv_tau;
		}
	}

	Lbi nc, nf;
	if (number_line > 0)
	{
		for (Lbi i=0; i<number_line; i++)
		{
			nc = NeborCF[i].nc;
			nf = NeborCF[i].nf;

			Line_Inf_C[i].areac_n = nc;
			Line_Inf_C[i].areaf_n = nf;

			Line_Inf_F[i].areac_n = nc;
			Line_Inf_F[i].areaf_n = nf;

			Line_init(NeborCF[i].direction, AreaDefine[nc], Areas[nc], Line_Inf_C[i], AreaDefine[nf], Areas[nf], Line_Inf_F[i]);
		}
	}

	if (number_line2 > 0)
	{
		for (Lbi i=0; i<number_line2; )
		{
			nc = NeborCC[i].nc;
			nf = NeborCC[i].nf;

			Line_Inf_CC[i].areac_n = nc;
			Line_Inf_CC[i].areaf_n = nf;

			Line_init(NeborCC[i].direction, AreaDefine[nc], Areas[nc], Line_Inf_CC[i], AreaDefine[nf], Areas[nf], Line_Inf_CC[i+1]);

			Line_Inf_CC[i+1].areac_n = NeborCC[i+1].nc;
			Line_Inf_CC[i+1].areaf_n = NeborCC[i+1].nf;

			P4 Line_tmp;
			Line_tmp = Line_Inf_CC[i+1].line_C;
			Line_Inf_CC[i+1].line_C = Line_Inf_CC[i+1].line_F;
			Line_Inf_CC[i+1].line_F = Line_tmp;

			//CC对应点数相同

			i = i+2;
		}
	}
}