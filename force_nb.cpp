//Program to calculate forces from Verlet neighbor list
double forcenb_vlst_u(double bl, int np, double **r, int npa,
		unsigned short int **vl, unsigned short int *nl,
		double *rc2, double *sig2, double eps)
{
	int i,j,k,lm,l,csize=1;
	double rx,ry,rz,ds2,ds2s,dsi6,
				 fijx,fijy,fijz,fij,pe;
	static double ***g; //force collectors and cut-off for each thread
	static double **sigma,**rcut2;
	static int m=0;

	//Allocate private space for each thread
	lm=omp_get_max_threads();

	if(m==0){
		sigma = new double*[lm];
		rcut2 = new double*[lm];
		g = new double**[lm];
		for(i=0;i<lm;i++){
			sigma[i] = new double[3];
			rcut2[i] = new double[3];
			g[i] = new double*[np];
			for(j=0;j<np;j++) g[i][j] = new double[4];
		}

		// Copy all the parameters to each thread
		for(i=0;i<lm;i++){
			for(k=0;k<3;k++){
				sigma[i][k]=sig2[k];
				rcut2[i][k]=rc2[k];
			}
		}	
		m++;
	}

	//Private read data for each thread
	for(i=0;i<lm;i++){
		for(j=0;j<np;j++){
			for(k=0;k<4;k++){
				if(k<2) g[i][j][k]=r[j][k];
				else g[i][j][k]=0.0;
			}
		}
	}

	/*
	 * bl - box length
	 * rx,ry,rz - distances
	 * ds - distance
	 * ds2 - sqare of the distance
	 * rc - square of potential cut off 
	 * pe - potential energy 
	 * vij - pair potential
	 * fij - pair force
	 * fijx, fijy, fijz - force components
	 */

 pe=0.0;  // Initialization of potential energy
 /*************OpenMP directive over the number of particles *************/
#pragma omp parallel for firstprivate(bl,np,npa,eps)\
	private(k,j,l,rx,ry,ds2,ds2s,dsi6,fij,fijx,fijy)\
	schedule(dynamic,csize)	reduction(+:pe)
	for(i=0;i<np-1;i++){

		l=omp_get_thread_num();

		for(k=0;k<nl[i];k++){

			j=vl[i][k];

			rx=g[l][j][0]-g[l][i][0];
			ry=g[l][j][1]-g[l][i][1];

			rx-=bl*rint(rx/bl);			
			ry-=bl*rint(ry/bl);

			ds2=rx*rx+ry*ry;

			if(i<npa && j<npa){
				if(ds2<=rcut2[l][0]){

					ds2s=sigma[l][0]/ds2;
					dsi6=ds2s*ds2s*ds2s;
					pe+=(4.0*eps*(dsi6-1.0)*dsi6)+eps;
					fij=24.0*eps*(2.0*dsi6-1.00)*dsi6/ds2;

					fijx=fij*rx;
					fijy=fij*ry;

					g[l][i][2]-=fijx;
					g[l][i][3]-=fijy;

					g[l][j][2]+=fijx;
					g[l][j][3]+=fijy;

				}
			}else if(i>=npa && j>=npa){
				if(ds2<=rcut2[l][1]){

					ds2s=sigma[l][1]/ds2;
					dsi6=ds2s*ds2s*ds2s;
					pe+=(4.0*eps*(dsi6-1.0)*dsi6)+eps;
					fij=24.0*eps*(2.0*dsi6-1.00)*dsi6/ds2;

					fijx=fij*rx;
					fijy=fij*ry;

					g[l][i][2]-=fijx;
					g[l][i][3]-=fijy;

					g[l][j][2]+=fijx;
					g[l][j][3]+=fijy;

				}
			}else{
				if(ds2<=rcut2[l][2]){

					ds2s=sigma[l][2]/ds2;
					dsi6=ds2s*ds2s*ds2s;
					pe+=(4.0*eps*(dsi6-1.0)*dsi6)+eps;
					fij=24.0*eps*(2.0*dsi6-1.00)*dsi6/ds2;

					fijx=fij*rx;
					fijy=fij*ry;

					g[l][i][2]-=fijx;
					g[l][i][3]-=fijy;

					g[l][j][2]+=fijx;
					g[l][j][3]+=fijy;

				}/* End of force cut-off */
			}/* End of particle identification IF */
		}/* End of inner loop */
	}/* End of outer loop */

	// Collect forces from all threads and assign it to the older array 
	for(i=0;i<np;i++)for(j=0;j<lm;j++)for(k=2;k<4;k++) r[i][k+2]+=g[j][i][k];

	return pe;
}
