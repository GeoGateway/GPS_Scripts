#include <math.h>
#define NRANSI
#include "nr_util.h"

float NR_elle(float phi, float ak)
{
	float rd(float x, float y, float z);
	float rf(float x, float y, float z);
	float cc,q,s;

	s=sin(phi);
	cc=NR_SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(NR_SQR(s*ak))*rd(cc,q,1.0)/3.0);
}
#undef NRANSI
