float NR_erff(float x)
{
	float NR_gammp(float a, float x);

	return x < 0.0 ? -NR_gammp(0.5,x*x) : NR_gammp(0.5,x*x);
}
