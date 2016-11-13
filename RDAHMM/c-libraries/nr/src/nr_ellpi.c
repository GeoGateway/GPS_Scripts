#include <math.h>
#define NRANSI
#include "nr_util.h"

float NR_ellpi(float phi, float en, float ak)
{
	float rf(float x, float y, float z);
	float NR_rj(float x, float y, float z, float p);
	float cc,enss,q,s;

	s=sin(phi);
	enss=en*s*s;
	cc=NR_SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-enss*NR_rj(cc,q,1.0,1.0+enss)/3.0);
}
#undef NRANSI
