/* Velocity-verlet integrator Langevin bath (Langevin Equation) */
//
double move(int np, double **r, double **zet, double **tht, double dt,
		double te, double bl, int npa, double rcut, double rcld2, double *rcut2,
		double *sig2, double eps, long int *seed, long int *stp, double *af,
		double *ptaf, int *npcrd1, int *nnb, int *api, int *nap, double gama)
{
	/*
	 * nn - positions
	 * crd - no. Coordinate of each particle
	 * np - number of particles
	 * i,j  - indices
	 */
	int i,j;
	static int k=0;
	double pe;
	static double dtb2,dt2b2,pf,dtsq,
								dt3b2,sigma,**ct;

	/* position r[..][0-1], velocity r[..][2-3], force r[..][4-5], 
	 * old force r[..][6-7]
	 * rad - sqrt(2dt)
	 * dt - time step 
	 * te - temperature
	 */

	if(k==0){
		pf=1.0/2.0/sqrt(3.0);
		dtsq=sqrt(dt);
		dtb2 = dt/2.0;
		dt2b2 = dt*dtb2;
		dt3b2=pow(dt,3.0/2.0);
		sigma=sqrt(2.0*te*gama); // k_B, friction coefficient and mass is unity
		ct = new double*[np];
		for(i=0;i<np;i++) ct[i]=new double[2];
		k++;
	}

	for(i=0;i<np;i++){

		ct[i][0]=dt2b2*(r[i][4]-gama*r[i][2])+sigma*dt3b2*(0.50*zet[i][0]+pf*tht[i][0]);
		ct[i][1]=dt2b2*(r[i][5]-gama*r[i][3])+sigma*dt3b2*(0.50*zet[i][1]+pf*tht[i][1]);

		r[i][0] = r[i][0]+r[i][2]*dt+ct[i][0];
		r[i][1] = r[i][1]+r[i][3]*dt+ct[i][1];
		
		// Assign forces at updated step to the previous step 
		r[i][6]=r[i][4];
		r[i][7]=r[i][5];
	}

	// Calculate forces using Verlet neighbour list at updated positions 
	pe=force(bl,np,r,npa,rcut,rcut2,sig2,eps);

	// Apply active force in random directions
	actforce(npa,r,af,ptaf,bl,dt,seed,stp,
			npcrd1,rcld2,nnb,api,nap);
	// if(*stp%500==0) cout << "No of active B particles = "<<" "<< nac << endl;

	for(i=0;i<np;i++){
		r[i][2] = r[i][2]+dtb2*(r[i][4]+r[i][6])-dt*gama*r[i][2]+sigma*dtsq*zet[i][0]-gama*ct[i][0]; 
		r[i][3] = r[i][3]+dtb2*(r[i][5]+r[i][7])-dt*gama*r[i][3]+sigma*dtsq*zet[i][1]-gama*ct[i][1]; 
	}

	return pe/(double)np;
}
