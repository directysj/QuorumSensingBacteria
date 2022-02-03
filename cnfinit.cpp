/* Generate initial coordinates for the 50:50 mixture */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

#include "ran3.cpp"
int main()
{
	long int idum;
	int i,j,k,l,nps,npl,np,x,is;
	double bl,sigaa,sigbb,sigab,de;

	fstream fp;

	// Configuration parameters
	nps=500,npl=500;
	de=0.40;
	sigaa=1.40,sigbb=1.0;
	sigab=0.50*(sigaa+sigbb);
	np=nps+npl;

	// One of the side of square (In case of 2D)
	bl=pow((double)np/de,1.0/2.0);
	//bl*=sigaa;

	cout <<"np = " << np <<"  " << "bl = " << bl <<"  " << "de = " << de << endl;

	/* Generate seed for the rand function */
	fp.open("seed.dat",ios::in);
	if(fp){
		fp >> is;
		fp.close();
	}else{
		cout << "ERROR: seed.dat file does not exist" << endl;
		exit(EXIT_FAILURE);
	}

	for(i=0;i<is;i++) x=rand();
	fp.open("seed.dat",ios::out);
	fp << x/100000<< endl;
	fp.close();
	idum=x;
	
	// cout << "SEED = " << idum << endl;
	// exit(EXIT_SUCCESS);

	int *idp,flg,id1,id2;
	double **r,dx,dy,ds2,ds,rmin;

	idp = new int[np];
	r = new double*[np];
	for(i=0;i<np;i++) r[i] = new double[2];

	for(i=0;i<np;i++){
		if(i<nps) idp[i]=1;
		else idp[i]=2;
	}

	//for(i=0;i<np;i++) cout << " idp[" << i <<"]="<< idp[i] << endl; 
	//cout << "CODE is running" << endl;
	//exit(EXIT_SUCCESS);

	flg = 0;
	for(i=0;i<np;i++){

		if(flg==0){

			k=0;
			// Generate random number between -1 and +1  
			r[i][0] = (2.0*ran3(&idum)-1.0)*bl;
			r[i][1] = (2.0*ran3(&idum)-1.0)*bl;

			flg=1;

		}else{

			do{

				id1=idp[i];
				r[i][0] = (2.0*ran3(&idum)-1.0)*bl;
				r[i][1] = (2.0*ran3(&idum)-1.0)*bl;

				for(l=k;l>=0;l--){

					id2=idp[l];
					dx = r[i][0] - r[l][0];
					dy = r[i][1] - r[l][1];

					/* Minimum image convention */
					dx = dx - bl*rint(dx/bl);
					dy = dy - bl*rint(dy/bl);

					ds2 = dx*dx + dy*dy;
					ds=sqrt(ds2);

					if(id1==1 && id2==1) rmin=1.40;  // AA
					else if (id1==2 && id2==2) rmin=1.00;  // BB
					else rmin=1.20;  // AB

					if(ds<rmin) break;

				}

			}while(ds<rmin);

			k++;

		}

			cout <<"Particle index = "<< i << endl;
		/* 
		 * Minimum image convention: to bring particles 
		 * in the original simulation box
		 */
		r[i][0] -= bl*rint(r[i][0]/bl); 
		r[i][1] -= bl*rint(r[i][1]/bl); 

		// cout << "Particle index = "<< i<< endl);

	}

	fp.open("icrd0p4.dat",ios::out);
	for(i=0;i<np;i++) fp << r[i][0] << "  "<< r[i][1] << endl;
	fp.close();

	return 0;
}
