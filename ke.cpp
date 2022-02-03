double kinen(int np, double **r)
{
	int i,j;
	double ke;

	ke=0.0;
	for(i=0;i<np;i++) ke+=(r[i][2]*r[i][2]+r[i][3]*r[i][3]);
	ke*=0.50;

	return ke;
}
