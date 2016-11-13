float NR_gammp(float a, float x)
{
	void NR_gcf(float *gammcf, float a, float x, float *gln);
	void NR_gser(float *gamser, float a, float x, float *gln);
	void NR_error(char error_text[]);
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) NR_error("Invalid arguments in routine NR_gammp");
	if (x < (a+1.0)) {
		NR_gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		NR_gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
