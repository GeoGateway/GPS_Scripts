#include <math.h>
#define EPS 1.0e-5
#define JMAX 20

float NR_qtrap(float (*func)(float), float a, float b)
{
	float NR_trapzd(float (*func)(float), float a, float b, int n);
	void NR_error(char error_text[]);
	int j;
	float s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=NR_trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		if (s == 0.0 && olds == 0.0 && j > 6) return s;
		olds=s;
	}
	NR_error("Too many steps in routine NR_qtrap");
	return 0.0;
}
#undef EPS
#undef JMAX
