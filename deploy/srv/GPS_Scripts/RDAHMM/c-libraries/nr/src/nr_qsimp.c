#include <math.h>
#define EPS 1.0e-6
#define JMAX 20

float NR_qsimp(float (*func)(float), float a, float b)
{
	float NR_trapzd(float (*func)(float), float a, float b, int n);
	void NR_error(char error_text[]);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=NR_trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		if (s == 0.0 && os == 0.0 && j > 6) return s;
		os=s;
		ost=st;
	}
	NR_error("Too many steps in routine NR_qsimp");
	return 0.0;
}
#undef EPS
#undef JMAX
