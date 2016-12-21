#include <math.h>
#define FACTOR 1.6
#define NTRY 50

int NR_zbrac(float (*func)(float), float *x1, float *x2)
{
	void NR_error(char error_text[]);
	int j;
	float f1,f2;

	if (*x1 == *x2) NR_error("Bad initial range in NR_zbrac");
	f1=(*func)(*x1);
	f2=(*func)(*x2);
	for (j=1;j<=NTRY;j++) {
		if (f1*f2 < 0.0) return 1;
		if (fabs(f1) < fabs(f2))
			f1=(*func)(*x1 += FACTOR*(*x1-*x2));
		else
			f2=(*func)(*x2 += FACTOR*(*x2-*x1));
	}
	return 0;
}
#undef FACTOR
#undef NTRY
