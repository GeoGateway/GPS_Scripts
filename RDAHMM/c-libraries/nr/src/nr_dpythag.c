#include <math.h>
#define NRANSI
#include "nr_util.h"

float NR_dpythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+NR_DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_DSQR(absa/absb)));
}
#undef NRANSI
