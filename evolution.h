
void evolution(Lbi n, Lbi *asp, Lbi *ls3cffcp, Lbi *ltcffcp, Lbi *lccp, Lbi level, dim3 block, dim3 *grid, AREA *areas, D2Q9 *pps_h, D2Q9 *pps_d, Inf_Ff_store *line_Inf_c, 
	Inf_Ff_store *line_Inf_f, Inf_Ff_store *line_Inf_cc, Ff_Store *line_t3, FM_spline3 *line_S3, cusparseHandle_t handle)
{
	for (Lbi i=0; i<mi[level]; i++)
	{
		if (level == LEVEL)
			return;
		if (!(level == 0 || i == 0))
		{	
			for (Lbi b=ltcffcp[level]; b<ltcffcp[level+1]; b++)
			{
				F_k3 *ft = &line_t3[b].f_t1;
				temporal_interpolation<<< int (ceil(line_Inf_c[b].Nf/double (_BLOCK_1D))), _BLOCK_1D>>>
					(areas[line_Inf_c[b].areaf_n], pps_d[line_Inf_c[b].areaf_n], line_Inf_c[b], 
					line_t3[b].N, ft[line_t3[b].marker[0]], ft[line_t3[b].marker[1]], ft[line_t3[b].marker[2]], (1.0*i/mi[level]));
			}
		}
		//output(1122, areas, pps_h, pps_d);
		for (Lbi a=asp[level]; a<asp[level+1]; a++)
		{
			evolution_area_stb( block, grid[a], areas[a], pps_d[a]);
		}

		if (number_line2 != 0)
		{
			for (Lbi c=lccp[level]; c<lccp[level+1]; c++)
				c_c_boundary_transfer<<< int (ceil(line_Inf_cc[c].Nc/double (_BLOCK_1D))), _BLOCK_1D>>>
				(areas[line_Inf_cc[c].areac_n], areas[line_Inf_cc[c].areaf_n], 
				pps_d[line_Inf_cc[c].areac_n], pps_d[line_Inf_cc[c].areaf_n], 
				line_Inf_cc[c]);
		}

		if (level != LEVEL-1)
		{
			for (Lbi b=ls3cffcp[level]; b<ls3cffcp[level+1]; b++)
			{
				space_interpolation_Spline3(areas[line_Inf_c[b].areac_n], areas[line_Inf_c[b].areaf_n], 
					pps_d[line_Inf_c[b].areac_n], pps_d[line_Inf_c[b].areaf_n], 
					line_Inf_c[b], line_t3[b], line_S3[b], handle);
				for (Lbi p=0; p<3; p++)
				{
					Lbi bb = line_t3[b].marker[p]-1; 
					if (bb == -1) line_t3[b].marker[p] = 2;
					else line_t3[b].marker[p] = bb;
				}
			}
		}
		evolution(n, asp, ls3cffcp, ltcffcp, lccp, level+1, block, grid, areas, pps_h, pps_d, line_Inf_c, line_Inf_f, line_Inf_cc, line_t3, line_S3, handle);
		if (level != 0 && i == mi[level]-1) 
		{
			for (Lbi b=ltcffcp[level]; b<ltcffcp[level+1]; b++)
			{
				F_k3 *ft = &line_t3[b].f_t1;
				c_f_boundary_transfer<<< int (ceil(line_Inf_c[b].Nf/double (_BLOCK_1D))), _BLOCK_1D>>>
					(areas[line_Inf_c[b].areaf_n], pps_d[line_Inf_c[b].areaf_n], 
					line_Inf_c[b], line_t3[b].N, ft[line_t3[b].marker[2]]);
				f_c_boundary_transfer<<< int (ceil(line_Inf_f[b].Nc/double (_BLOCK_1D))), _BLOCK_1D>>>
					(areas[line_Inf_f[b].areac_n], areas[line_Inf_f[b].areaf_n], 
					pps_d[line_Inf_f[b].areac_n], pps_d[line_Inf_f[b].areaf_n], 
					line_Inf_f[b]);
			}
		}
	}
}