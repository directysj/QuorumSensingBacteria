/*
 * Radial distribution function
 */
void rdf(double bl, int np, double de, 
		double **r, int nap, int *api)
{
	/*
	 * i,j,k,jj indices
	 * mbil - max. bin length
	 * np - number of particles
	 */

	int i,j,bn,tn,l,mt,m,n;
	static int mbil,csize,npl,nps;
	static unsigned long int k=0,p=0;
	static int **apig;

	mt=omp_get_max_threads();
	/*
	 * cvn - nor. constant
	 * bl - box length
	 * dr - size of a bin
	 * dis1,dis2,dis3,dis4 - temp var.s
	 * cbn - norm. for the bins
	 * gr - rad. dis. fn.
	 * ds2 - ds*ds
	 * ds  - distance
	 */

	double dis1,dis2,dis3,dis4,dis5,dis6;
	double ds2,ds,dx,dy,grac;
	static double ***g,*cbn,**gr,***lgr,
		      dr=0.01,pi,A;

	// Normalization constants for bins
	if(k==0){
		A = bl*bl;
		npl = np/2; // no of larger particles (50:50 mixture)
		nps = npl;
		csize=1;
		pi=4.0*atan(1.0);
		mbil=(int)(bl/2.0/dr);
		cbn = new double[mbil];
		gr = new double*[mbil];
		for(i=0;i<mbil;i++) gr[i] = new double[5];
		lgr = new double**[mt];
		apig = new int*[mt];
		g = new double**[mt];
		for(i=0;i<mt;i++){
			apig[i]=new int[npl];
			lgr[i]=new double*[mbil];
			g[i]=new double*[np];
			for(j=0;j<np;j++) g[i][j]=new double[2];
			for(j=0;j<mbil;j++) lgr[i][j]=new double[5];
		}

		for(i=0;i<mbil;i++){
			cbn[i]=0.0;
			for(j=0;j<5;j++) gr[i][j] = 0.0;
		}

		dis1=dr;
		dis2=0.0;
		for(j=0;j<mbil;j++){
			cbn[j]=pi*(dis1*dis1-dis2*dis2);
			dis1+=dr;
			dis2+=dr;
		}
	}

	k++;
	if(nap>=2) p++;
	// Allocate co-ordinates to all the threads
	for(j=0;j<mt;j++){
		for(i=0;i<nap;i++) apig[j][i]=api[i];
		for(i=0;i<np;i++)for(l=0;l<2;l++) g[j][i][l]=r[i][l];
		for(i=0;i<mbil;i++)for(l=0;l<5;l++) lgr[j][i][l]=0.0;
	}

#pragma omp parallel for firstprivate(np,nap,dr,mbil,npl)\
	private(j,dx,dy,ds2,ds,bn,tn,m,n) schedule(dynamic,csize)
	for(i=0;i<np;i++){

		tn=omp_get_thread_num();

		for(j=0;j<np;j++){

			if(nap>=2){
				if(i<nap && j<nap){

					m=apig[tn][i];
					n=apig[tn][j];

					if(m!=n){

						dx=g[tn][m][0]-g[tn][n][0];
						dy=g[tn][m][1]-g[tn][n][1];

						dx-=bl*rint(dx/bl);			
						dy-=bl*rint(dy/bl);			

						ds2=dx*dx+dy*dy;
						ds=sqrt(ds2);

						if(ds<0.25){
							cout << "grs.c:act "<<" "<< ds <<" "<< dx <<" "<< dy << endl;
							cout << "grs.c:act "<<" "<< m <<" "<< n << endl;
							exit(EXIT_FAILURE);
						}

						bn=(int)(ds/dr);
						if(bn<mbil) lgr[tn][bn][4] += 2.0;
					}
				}
			}

			if(i!=j){

				dx=g[tn][i][0]-g[tn][j][0];
				dy=g[tn][i][1]-g[tn][j][1];

				dx-=bl*rint(dx/bl);			
				dy-=bl*rint(dy/bl);			

				ds2=dx*dx+dy*dy;
				ds=sqrt(ds2);

				if(ds<0.25){
					cout << "grs.c:act "<< ds <<" "<< dx <<" "<< dy << endl;
					cout << "grs.c:act "<< m <<" "<< n << endl;
					exit(EXIT_FAILURE);
				}

				bn=(int)(ds/dr);

				if(bn<mbil){
					if(i<npl && j<npl) lgr[tn][bn][0] += 1.0;
					if(i>=npl && j>=npl) lgr[tn][bn][1] += 1.0;
					if(i<npl && j>=npl) lgr[tn][bn][2] += 1.0;
					if(i>=npl && j<npl) lgr[tn][bn][3] += 1.0;
				}
			}
		}
	}

	// Collect data from all threads
	for(i=0;i<mbil;i++){
		grac=0.0;
		for(j=0;j<mt;j++){
			for(l=0;l<4;l++) gr[i][l]+=lgr[j][i][l];
			if(nap>=2) grac += lgr[j][i][4];
		}
		if(nap>=2){
			grac /= ((double)(nap*(nap-1)));
			gr[i][4] += grac;
		}
	}

	ofstream out;

	if(k%100==0){
		out.open("rdf.dat",ios::out);
		for(i=0;i<mbil;i++){
			dis1=((double)i+0.5)*dr;
			dis2=(A*gr[i][0])/cbn[i]/((double)(k*npl*(npl-1)));
			dis3=(A*gr[i][1])/cbn[i]/((double)(k*nps*(nps-1)));
			dis4=(A*gr[i][2])/cbn[i]/((double)(k*(npl-1)*nps));
			dis5=(A*gr[i][3])/cbn[i]/((double)(k*(nps-1)*npl));
			dis6=(A*gr[i][4])/cbn[i]/(double)p;
			out << dis1 <<" "<< dis2 <<" "<< dis3 <<" "<< dis4 <<" "<< dis5 <<" "<< dis6 <<" "<< endl;
		}
		out.close();
	}

}
