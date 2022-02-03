/* Give random displacements according to Maxwell-Boltzmann distributions */
void grand_disp(int np, double **r, int dim, long int *iseed)
{
	int i,j;
	static int k=0;
	static long int idum;
	int is,x;

	for(i=0;i<np;i++)for(j=0;j<dim;j++) r[i][j]=gaus(iseed);

}
