
void output(Lbi n, AREA *areas, D2Q9 *pps_h, D2Q9 *pps_d)
{
	for (Lbi a=0; a<number_area; a++)
	{
		Lbui N = areas[a].nx*areas[a].ny;
		cudaMemcpy(pps_h[a].solid, pps_d[a].solid, N*sizeof(Lbi), cudaMemcpyDeviceToHost);
		Lbd* *Vh = &pps_h[a].x;
		Lbd* *Vd = &pps_d[a].x;
		for (Lbi b=0; b<25; b++)
		{
			cudaMemcpy(Vh[b], Vd[b], N*sizeof(Lbd), cudaMemcpyDeviceToHost);
		}
	}

	int ixy = 0;
	FILE *steady = NULL;
	char filename[50];
	sprintf(filename,"steady_%d.dat",n);
	steady=fopen(filename,"w");
	fprintf(steady,"Title= \"LBM 2D\"\n");
	fprintf(steady,"VARIABLES= \"X\",\"Y\",\"solid\",\"rho\",\"U\",\"V\",\n");

	for (Lbi a=0; a<number_area; a++)
	{
		fprintf(steady,"ZONE T= \"zone_level%d\",I= %d,J= %d,F=POINT\n",a, areas[a].nx, areas[a].ny);	
		for(int j=0;j<areas[a].ny;j++){
			for(int i=0;i<areas[a].nx;i++){
				ixy =  i+j*areas[a].nx;
				//if (pps_h[a].solid[ixy] > 100) system("pause");
				fprintf(steady,"%g\t%g\t%d\t%g\t%g\t%g\n",
					pps_h[a].x[ixy], pps_h[a].y[ixy], pps_h[a].solid[ixy], 
					pps_h[a].rho[ixy], pps_h[a].ux[ixy], pps_h[a].uy[ixy]);	
			}
		}
	}
	fclose(steady);
}