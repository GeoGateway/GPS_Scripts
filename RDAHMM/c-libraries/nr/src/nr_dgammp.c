double NR_dgammp(double a, double x)
{
	void NR_dgcf(double *gammcf, double a, double x, double *gln);
	void NR_dgser(double *gamser, double a, double x, double *gln);
	void NR_error(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) NR_error("Invalid arguments in routine NR_gammp");
	if (x < (a+1.0)) {
		NR_dgser(&gamser,a,x,&gln);
		return gamser;
	} else {
		NR_dgcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
