/* Program to generate Gaussian random deviates*/
double gaus(long int *idum)
{
	int i,j;
	long int seed;
	seed=*idum;
	double a1,a3,a5,a7,a9;
	a1=3.949846138,a3=0.252408784,a5=0.076542912;
	a7=0.008355968,a9=0.029899776;

	double sum,R,Rsq,gauss;
	sum=0.0;
	
	for(i=0;i<12;i++) sum += ran3(&seed);

	R = (sum-6.0)/4.0;
	Rsq = R*R;
	gauss = ((((a9*Rsq+a7)*Rsq+a5)*Rsq + a3)*Rsq + a1)*R;

	return gauss;
}
