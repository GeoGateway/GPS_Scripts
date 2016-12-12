#include <math.h>
#define NRANSI
#include "nr_util.h"

float NR_ellf(float phi, float ak)
{
	float rf(float x, float y, float z);
	float s;

	s=sin(phi);
	return s*rf(NR_SQR(cos(phi)),(1.0-s*ak)*(1.0+s*ak),1.0);
}
#undef NRANSI
