#include <math.h>

float NR_beta(float z, float w)
{
	float NR_gammln(float xx);

	return exp(NR_gammln(z)+NR_gammln(w)-NR_gammln(z+w));
}
