#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void NR_gcf(float *gammcf, float a, float x, float *gln)
{
	float NR_gammln(float xx);
	void NR_error(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=NR_gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) NR_error("a too large, ITMAX too small in NR_gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
