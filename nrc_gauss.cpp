/*
 * Gaussian random deviates from the Book 
 * "Numerical Recipies in C".
 */

void grand(int np, int crd, double **r)
{
	static int iset=0,i,j,ctr=0;
	static long idum;
	int x,is;
	static double gset;
	double fac,rsq,v1,v2;
	fstream fp;

	if(ctr==0){
		fp.open("seed.dat",ios::in);
		fp >> is;
		fp.close();
		for(i=0;i<is;i++)x=rand();

		fp.open("seed.dat",ios::out);
		fp << x/1000000 << endl;
		fp.close();

		ctr++;
		idum=x;
	}

	for(i=0;i<np;i++){
		for(j=12;j<15;j++){  // This index should be updated according to purpose

			if (idum < 0) iset=0; 
			if (iset == 0) { 
				do {
					v1=2.0*ran3(&idum)-1.0; 
					v2=2.0*ran3(&idum)-1.0; 
					rsq=v1*v1+v2*v2;
				} while (rsq >= 1.0 || rsq == 0.0); 

				fac=sqrt(-2.0*log(rsq)/rsq);
				gset=v1*fac;

				iset=1; 

				r[i][j]=v2*fac;

			} else { 
				iset=0; 
				r[i][j]=gset; 
			}
		}
	}
}
