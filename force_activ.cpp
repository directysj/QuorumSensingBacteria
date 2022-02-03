/* Function to introduce activity in the system */
void actforce(int npa, double **r, double *af, double *ptaf, double bl,
		double dt, long int *iseed, long int *st, int *npcrd1, double rcld2, 
		int *nnb, int *api, int *nap)
{
	/*
	 * npcrd1 - no of particles in the first coordination shell
	 */
	static int k=0,ncrd,ntaup,npb,np;
	int i,j,l=0;
	static double **ranforc,f0,taup;
	double ranx,rany,ranm;

	if(k==0){
		// As a 50:50 mixture npa and npb are equal
		npb = npa;
		np = npa + npb;
		f0 = *af; // Magnitude of active force
		taup = *ptaf; // Persistence time for active force
		ntaup = (int)(taup/dt);
		ranforc = new double*[npb];
		for(i=0;i<npb;i++) ranforc[i] = new double[2];
		k++;
	}

	if((*st)%ntaup==0){
		for(i=0;i<npb;i++) nnb[i]=0;
		neborlist(npa,np,r,bl,rcld2,nnb);
	}

	for(i=npa;i<np;i++){

		if(nnb[i-npa] >= *npcrd1){

			if((*st)%ntaup==0){
				// Assign the random force direction at each
				// persistence time for the active particles 
				ranx = 2.0*ran3(iseed)-1.0;
				rany = 2.0*ran3(iseed)-1.0;
				
        ranm = sqrt(ranx*ranx+rany*rany);

				ranforc[i-npa][0] = ranx/ranm;
				ranforc[i-npa][1] = rany/ranm;

				// Store the ID's of active particles 
				api[i-npa] = i;
			}

			r[i][4] += f0*ranforc[i-npa][0];
			r[i][5] += f0*ranforc[i-npa][1];

			l++;

		}

	} /* End of active particles loop */

	*nap = l;
	
}
