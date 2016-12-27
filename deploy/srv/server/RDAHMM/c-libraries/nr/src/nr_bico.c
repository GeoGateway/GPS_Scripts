#include <math.h>

float NR_bico(int n, int k)
{
	float NR_factln(int n);

	return floor(0.5+exp(NR_factln(n)-NR_factln(k)-NR_factln(n-k)));
}
