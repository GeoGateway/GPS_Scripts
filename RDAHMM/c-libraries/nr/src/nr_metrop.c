#include <math.h>

int NR_metrop(float de, float t)
{
	float NR_ran3(long *idum);
	static long gljdum=1;

	return de < 0.0 || NR_ran3(&gljdum) < exp(-de/t);
}
