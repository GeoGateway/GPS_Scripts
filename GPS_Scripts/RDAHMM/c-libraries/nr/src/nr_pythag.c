#include <math.h>
#define NRANSI
#include "nr_util.h"

float NR_pythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+NR_SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_SQR(absa/absb)));
}
#undef NRANSI
