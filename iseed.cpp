long int seedrv()
{
	int i;
	long int is,x,idum;

	/* Generate seed of the rand */
	fstream fp;

	fp.open("seed.dat",ios::in);

	if(fp){
		fp >> is ;
		fp.close();
	}else{
		cout << "iseed.cpp: ERROR in file seed.dat reading" << endl;
		exit(EXIT_FAILURE);
	}

	for(i=0;i<is;i++) x=rand();

	fp.open("seed.dat",ios::out);
	if(fp){
		fp << x/100000 << endl;
		fp.close();
	}else{
		cout << "iseed.cpp: ERROR in file seed.dat writing" << endl;
		exit(EXIT_FAILURE);
	}

	idum=x;

	return idum;
}
