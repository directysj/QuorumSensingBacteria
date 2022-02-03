/* Force with Verlet neighbor list */
#include "force_nb_vlst.cpp"
#include "force_nb.cpp"

#define max(a,b) (((a)>(b))?(a):(b)) 
//
double force(double bl, int np, double **r,
		 int npa, double rcut, double *rc2,
		 double *sig2, double eps)
{ 
	int i,j,k;
	static int k1=0;

	double dmax,a,pe;

	const double skn=1.0;
	static double **g,rlist;
	static unsigned short int **vl,*nl,updt=1;

	if(k1==0){
		nl = new unsigned short int[np];
		g = new double*[np];
		vl = new unsigned short int*[np];
		for(i=0;i<=np;i++){
			g[i] = new double[2];
			vl[i] = new unsigned short int[150];
		}

		k1++;
		rlist=rcut+1.0;
		rlist*=rlist;
	}

	/************ Initialization of force **************/

	for(i=0;i<np;i++)for(j=4;j<6;j++) r[i][j]=0.0;

	 // Check for the necessity of update

	dmax=0.0;

	for(i=0;i<np;i++){
		dmax=max(dmax,fabs(r[i][0]-g[i][0]));
		dmax=max(dmax,fabs(r[i][1]-g[i][1]));
	}

	dmax=dmax*3.4641;/* 2*sqrt(3*dmax^2) ??*/
	updt=skn<dmax;


	/********************* Main Force Loop ********************/

	if(updt){

		// Save current configuration
		for(i=0;i<np;i++)for(j=0;j<2;j++) g[i][j]=r[i][j];

		//Generate neighbor list
		//In future code this will be replaced by cell list 
		pe=forcenb_vlst(bl,np,r,npa,vl,nl,rc2,rlist,sig2,eps);

		//s1++;
	
	}else{

		//Use neighbor list to calculate 
		pe=forcenb_vlst_u(bl,np,r,npa,vl,nl,rc2,sig2,eps);

		//s2++;

	}/*End of non-bonding forces*/

	//k1++;
	//if(k1%1000)printf("%lf %lf %lf\n",s1/(s1+s2),s1,s2);
	
  return pe;
}/*End of force*/
