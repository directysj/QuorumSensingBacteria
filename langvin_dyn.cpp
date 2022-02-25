/* Langevin dynamics simulation of 50:50 Binary Mixture */
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <cstring>

using namespace std; 

#include "ran3.cpp"
#include "simparin.cpp"
//#include "intcnf.cpp"
#include "gauss.cpp"
#include "iseed.cpp"
#include "veldistb.cpp"
#include "nnb.cpp"
#include "force_activ.cpp"
#include "force.cpp"
#include "gausranlang.cpp"
#include "mv_vvlt_ld.cpp"
#include "ke.cpp"
#include "trjrst.cpp"
#include "grs.cpp"

int main()
{
	/*
	 * ns - number of steps
	 * tf - Final time of the simulation
	 * nsi - number of sampling interval
	 * si - sampling interval in time
	 * neq - equilibration steps
	 * teq - equilibration time
	 * np - number of particles
	 * dim - no of spatial dimensions
	 * rcld - Cutoff distance for the active particles
	 * nap - No of active particles
	 * api - Active particles' ID's
	 */ 
	int j,n,flg,nsi,np,npa,npb,neq,resi,
			nav,npcrd1,*nnb,dim,nap,*api;
	long int i,ns,nsum_runvar,seed,seed1;
			 
	/*
	 * te - temperature
	 * de - density
	 * dt - time step
	 * bl - box length
	 * v - P.E.
	 * d - diffusion constant
	 * rad - sqrt(2dt)
	 * rcut - P.E. Cut off
	 */
	double te,de,dt,bl,tf,si,teq,af,
				 siga,sigb,eps,pe,ke,etot,tinst; 
	double pesum,kesum,etotsum,tmpsum,ptaf,
				 pesqsum,kesqsum,etotsqsum,tmpsqsum;
	double peav,keav,etotav,tmpav,pesqav,kesqav,
				 etotsqav,tmpsqav,flpe,flke,fletot,
				 fltmp,napsum,napsqsum,napav,napsqav,
				 flnap,gama;

	/*
	 * nn - postions 
	 * crd - no. coordinate of each particle
	 * position r[..,0-1],velocity r[..,2-3],
	 * force r[..,4-5], old force r[..][6-7]
	 * two uncorrelated random distributions
	 * zet[i][0-1] and tht[i][0-1]
	 */
	double **r,tc,*sig2,*rcut2,
				 **zet,**tht,rcut,rcld,
				 rcld2;

	fstream inp;

	dim = 2;

	/* Read Coordinates */
	inp.open("inp.dat",ios::in);
	if(inp){
		inp >> npa >> npb >> nav >> sigb >> eps >> tf >> si >> teq >> te >> de >> dt >> gama;
	inp.close();
	}else{
		cout <<"inp.dat :: Input parameters file opening ERROR"<<endl;
		exit(EXIT_FAILURE);
	}

	// Read active force parameters
	inp.open("activ.dat",ios::in);
	if(inp){
	inp >> af >> ptaf >> rcld >> npcrd1;
	inp.close();
	}else{
		cout <<"active.dat :: Input parameters file opening ERROR"<<endl;
		exit(EXIT_FAILURE);
	}

	np=npa+npb;

	// cout << "np = " << np << " " << "Active force = " << af << endl;

	bl = pow((double)np/de,1.0/(double)dim);
	rcut = pow(2.0,1.0/6.0);
	rcld2 = rcld*rcld;

	nnb = new int[npb];
	api = new int[npb];
	r = new double*[np];
	zet = new double*[np];
	tht = new double*[np];
	for(i=0;i<np;i++){
		r[i] = new double[4*dim];
		zet[i] = new double[dim];
		tht[i] = new double[dim];
	}

	sig2 = new double[3];
	rcut2 = new double[3];

	// Initialize no of nearest neighbours and the indices of 
	// active particles
	for(i=0;i<npb;i++) nnb[i]=api[i]=0;

	inp.open("pcheck.dat",ios::out);
	inp << "Particles= "<<  np << endl;
	inp << "Final time= "<< tf << endl;
	inp << "Density= "<< de << endl;
	inp << "Temperature= "<< te << endl;
	inp << "Box Length= "<< bl << endl;
	inp.close();

	// Simulation parameters
	siga=simparinit(sigb,sig2,rcut,rcut2);

	//pos=fopen("icrd0p6.dat","r");
	inp.open("icrd.dat",ios::in);
	if(inp){
	for(i=0;i<np;i++) inp >> r[i][0] >> r[i][1] >> r[i][2] >> r[i][3];
	inp.close();
	}else{
		cout << "ERROR:: Initial configuration could not be read or does not exist"<< endl;
		exit(EXIT_FAILURE);
	}

	// Apply minimum image convention to ensure the particles 
	// in the central simulation box
	for(i=0;i<np;i++)for(j=0;j<dim;j++) r[i][j]-=bl*rint(r[i][j]/bl);

	// Generate initial configuration on a FCC Lattice
	//init(np,npa,siga,sigb,bl,de,r);

	// Generate two uncorrelated Gaussian random variables
	seed=seedrv();
	seed1=seedrv();

	// Random Gaussian distribution of velocities
	veldismb(np,r,te,dim,&seed);

	nsi=(int)(si/dt);  // No of sampling interval
	neq=(int)(teq/dt); // No of equilibration time
	tf+=teq;

	nap = 0; // Initial no of active particles
	nsum_runvar=0; // Running variable summation number
	tc=0.0; // Current time initialization
	i=0; // Step number initialization

	pesum=kesum=etotsum=tmpsum=napsum=0.0;
	pesqsum=kesqsum=etotsqsum=tmpsqsum=napsqsum=0.0;

	// Restart the code if required that is decided by 0/1
	traj(&tc,&i,np,npb,r,&nap,api,0,dt);

	/********************* Main Loop *********************/
	while(tc<=tf){

		if(i==0){
			// Compute forces using Verlet neighbour list
			pe=force(bl,np,r,npa,rcut,rcut2,sig2,eps);
			cout << "Initial potential energy per particle = " << pe/(double)np <<endl;
			cout << "Steps    Time    P.E.   K.E.   Energy   Temperature"<< endl;
		}

		/*********** Uncorrelated Gaussian random deviates for the displacement ***********/
		// Random numbers Zeta(t) for each particle
		grand_disp(np,zet,dim,&seed);

		// Random numbers theta(t) for each particle
		grand_disp(np,tht,dim,&seed1); 
		/*********** Uncorrelated Gaussian random deviates for the displacement ***********/

		// Move the particles using Velocity verlet with Langevin noise
		pe=move(np,r,zet,tht,dt,te,bl,npa,rcut,rcld2,rcut2,sig2,
				eps,&seed,&i,&af,&ptaf,&npcrd1,nnb,api,&nap,gama);

		// Calculate kinetic energy of the system
		ke=kinen(np,r);

		// Instantaneous temperature
		tinst = ke/(((double)np-1.0));

		// Per particle K.E.
		ke/=(double)np;

		// Total energy of the system
		etot=ke+pe;

		// if(i%5000==0) cout << i << " " << (double)i*dt << " " << pe<< " " << ke << " " << etot<< " " << tinst << endl;

		// Storing trajectory and calculation of correlations
		if(i>neq){ 

			if(i%nsi==0){  
				// store trajectories
				 traj(&tc,&i,np,npb,r,&nap,api,0,dt);

				// calculate radial distribution functions
				// rdf(bl,np,de,r,nap,api);
			}

			// Accumulating running variables
			pesum = pesum + pe;
			kesum = kesum + ke;
			etotsum = etotsum + etot;
			tmpsum = tmpsum + tinst;
			napsum = napsum + (double)nap;

			pesqsum = pesqsum + pe*pe;
			kesqsum = kesqsum + ke*ke;
			etotsqsum = etotsqsum + etot*etot;
			tmpsqsum = tmpsqsum + tinst*tinst;
			napsqsum = napsqsum + (double)(nap*nap);

			nsum_runvar+=1;

			if(i%nav == 0){

				peav = pesum/(double)nsum_runvar;
				keav = kesum/(double)nsum_runvar;
				etotav = etotsum/(double)nsum_runvar;
				tmpav = tmpsum/(double)nsum_runvar;
				napav = napsum/(double)nsum_runvar;

				pesqav = pesqsum/(double)nsum_runvar;
				kesqav = kesqsum/(double)nsum_runvar;
				etotsqav = etotsqsum/(double)nsum_runvar;
				tmpsqav = tmpsqsum/(double)nsum_runvar;
				napsqav = napsqsum/(double)nsum_runvar;

				flpe = sqrt(fabs(pesqav-peav*peav));
				flke = sqrt(fabs(kesqav-keav*keav));
				fletot = sqrt(fabs(etotsqav-etotav*etotav));
				fltmp = sqrt(fabs(tmpsqav-tmpav*tmpav));
				flnap = sqrt(fabs(napsqav-napav*napav));

				cout <<"After "<< i << " steps"<< endl;
				cout <<"Averaged over "<< nsum_runvar << "steps"<< endl;

				cout <<"<PE> = " << peav <<  "+/- "<< flpe << endl;
				cout <<"<KE> = " << keav <<  "+/- "<< flke << endl;
				cout <<"<Etot> = " << etotav << "+/- "<<fletot<< endl;
				cout <<"<T> = "<< tmpav << "+/- "<< fltmp<< endl;
				cout <<"<NAP> = "<< (int)napav << "+/- "<< (int)flnap<< endl;

			}
		}

		tc+=dt; // Update current time
		i++;  // Update step number

	} 

	// store trajectories for last step numbers those could not
	// be saved because of maximum no of records of a file
	 traj(&tc,&i,np,npb,r,&nap,api,1,dt);
  
	inp.open("fconf.dat",ios::out);
	for(i=0;i<np;i++) inp <<r[j][0]<<" "<<r[j][1]<<" "<<r[j][2]<<" "<<r[j][3]<<endl;
	inp.close();


	delete [] nnb;
	delete [] api;
	delete [] r;
	delete [] zet;
	delete [] tht;
	delete [] sig2;
	delete [] rcut2;
	
	return 0;
}
