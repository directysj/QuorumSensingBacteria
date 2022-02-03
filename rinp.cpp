//Code to read the l'th file in the trajectory
int readtrj(int scnf, int np, int *nap, int **api, 
		int l, long *stp, double *tc, double ***g)
{
	int x,j,m;
	static char fname[20],fno[5],pre[]="trj/t",post[]=".bdt";
	ifstream inp;

	//Generate file name of the next record

	sprintf(fno,"%d",l);
	strcpy(fname,pre);
	strcat(fname,fno);
	strcat(fname,post);
	cout << fname << endl;

	//Print out the next record name
	//Open the next trajectory file that is to be analyzed

	inp.open(fname,ios::binary | ios::in);

	if(inp){

		// Read a full trajectory file
		for(x=0;x<scnf;x++){
			// Steps
			inp.read((char*) &stp[x],sizeof(long int));
 
			// Current time 
			inp.read((char*) &tc[x], sizeof(double));
			
			//Positions of particles in a configuration
			for(j=0;j<np;j++) inp.read((char*) g[x][j],4*sizeof(double));  
			
			// Active particles' information: No and ID's
			inp.read((char*) &nap[x],sizeof(int));
			inp.read((char*) api[x],nap[x]*sizeof(int));

		}
		//Close the current trajectory file
		inp.close();

	}else{
		cout << "Error: FILE could not be opened or does not exist"<< endl;
		return 0;
	}

	//Finally return opened number of records before final exit
	return x;
}
