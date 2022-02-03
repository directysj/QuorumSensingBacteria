// Initial configuration on a HCP lattice in 2D
void init(int np, int npa, double siga, double sigb, 
		double bl, double de, double **r)
{
	int i,j,k,l,n;
	// np - no of particles in the system
	// l - for particle index
	// bl - length of the square
	// de - Density
	// psy - Planar separation in Y direction
	double psy;
	int lx,ly,lxn,lyn;
	FILE *fp;

	k=100000;

	psy = 0.50*pow(3.0,1.0/2.0);

	// Initial coordinates for the Larger particles
	lx=(int)(bl/siga);
	ly=lx;

	// cout <<"bl = " << bl <<" "<<"Lx = "<< lx <<" "<<" Ly = "<< ly<< endl;

	l=0;	
	for(j=0;j<ly;j++){
		for(i=0;i<lx;i++){
			if(j%2==0){
				r[l][0]=i*siga;
				r[l][1]=j*psy*siga;
			}else{
				r[l][0]=i*siga+0.50*siga;
				r[l][1]=j*psy*siga;
			}
			l++;
			if(l==npa){
				k=j;
				n=i;
				break;
			}
		}
		if(j==k) break;
	}

	// cout <<"k = "<< k<< endl;

	// Initial coordinates for the Smaller particles
	lxn=(int)(bl/sigb);
	lyn=lxn;

	// cout <<"New Lx =" << lxn << "New Ly = " << lyn <<endl);

	for(j=k+1;j<lyn;j++){
		for(i=0;i<lxn;i++){
			if(j%2==0){
				r[l][0]=i*sigb;
				r[l][1]=k*psy*siga+(j-k)*psy*sigb;
			}else{
				r[l][0]=i*sigb+0.50*sigb;
				r[l][1]=k*psy*siga+(j-k)*psy*sigb;
			}
			l++;
			if(l>=np) break;
		}
		if(l>=np) break;
	}

	fp.open("crdhcp.dat",ios::out);
	for(i=0;i<np;i++) fp << r[i][0] <<" "<< r[i][1]<< endl;
	fp.close();

}
