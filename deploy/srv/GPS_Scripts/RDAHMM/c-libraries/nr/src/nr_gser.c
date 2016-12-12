#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void NR_gser(float *gamser, float a, float x, float *gln)
{
	float NR_gammln(float xx);
	void NR_error(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=NR_gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) NR_error("x less than 0 in routine NR_gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		NR_error("a too large, ITMAX too small in routine NR_gser");
		return;
	}
}
#undef ITMAX
#undef EPS
