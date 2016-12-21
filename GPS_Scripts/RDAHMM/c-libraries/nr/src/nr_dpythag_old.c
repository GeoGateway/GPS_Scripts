#include <math.h>
#define NRANSI
#include "nr_util.h"

double NR_dpythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+NR_DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_DSQR(absa/absb)));
}
#undef NRANSI
