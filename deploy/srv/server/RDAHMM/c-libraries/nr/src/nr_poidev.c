#include <math.h>
#ifndef PI
#ifdef M_PI
#define PI  M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif

float NR_poidev(float xm, long *idum)
{
	float NR_gammln(float xx);
	float NR_ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= NR_ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-NR_gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*NR_ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-NR_gammln(em+1.0)-g);
		} while (NR_ran1(idum) > t);
	}
	return em;
}
#undef PI
