// List the particles those will be active
// according to the no of nearest neighbors
void neborlist(int npa, int np, double **r, 
		double bl, double rcld2, int *nnb)
{
	int i,j,tn,n;
	double ds2,dx,dy;
	static double ***g;
	static int k=0,mt,dim,**nbl;

	if(k==0){
		dim=2;
		mt=omp_get_max_threads();
		g = new double**[mt];
		nbl = new int*[mt];
		for(i=0;i<mt;i++){
			g[i] = new double*[np];
			nbl[i] = new int[npa];
			for(j=0;j<np;j++) g[i][j] = new double[dim];
		}
		k++;
	}

	// Copy data to each slave thread and initialization
	for(i=0;i<mt;i++){
		for(j=0;j<np;j++){
			if(j<npa) nbl[i][j] = 0;
			for(n=0;n<dim;n++) g[i][j][n] = r[j][n];
		}
	}

#pragma omp parallel for firstprivate(npa,np,bl,rcld2)\
	private(j,tn,dx,dy,ds2) schedule(dynamic,1)
	for(i=npa;i<np;i++){

		tn=omp_get_thread_num();
		for(j=0;j<np;j++){

			if(i!=j){

				dx = g[tn][i][0] - g[tn][j][0];
				dy = g[tn][i][1] - g[tn][j][1];

				dx-=bl*rint(dx/bl);
				dy-=bl*rint(dy/bl);

				ds2 = dx*dx + dy*dy;

				if(ds2<=rcld2) nbl[tn][i-npa]+=1;

			}
		}
	}

	// Collect data from all the threads and copy into the array
	for(i=npa;i<np;i++)for(j=0;j<mt;j++) nnb[i-npa] += nbl[j][i-npa];

}
