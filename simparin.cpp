double simparinit(double sigbb, double *sig2,
		 double rcut, double *rcut2)
{
	int i;

	double rcaa,rcbb,rcab,dsqi,dsi6,dsi12;
	double sigaa,sigab;

	// Particle diameters
	sigaa = 1.40*sigbb;  // sigma_AA
	sigab = 0.50*(sigaa+sigbb);  // sigma_AB

	sig2[0] = sigaa*sigaa;
	sig2[1] = sigbb*sigbb;
	sig2[2] = sigab*sigab;

	rcaa = rcut*sigaa;
	rcbb = rcut*sigbb;
	rcab = rcut*sigab;

	rcut2[0] = rcaa*rcaa;
	rcut2[1] = rcbb*rcbb;
	rcut2[2] = rcab*rcab;

	return sigaa;

}
