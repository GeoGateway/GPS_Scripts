float NR_erffc(float x)
{
	float NR_gammp(float a, float x);
	float NR_gammq(float a, float x);

	return x < 0.0 ? 1.0+NR_gammp(0.5,x*x) : NR_gammq(0.5,x*x);
}
