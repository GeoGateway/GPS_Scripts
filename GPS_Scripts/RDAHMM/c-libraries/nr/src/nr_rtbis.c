#include <math.h>
#define JMAX 40

float NR_rtbis(float (*func)(float), float x1, float x2, float xacc)
{
	void NR_error(char error_text[]);
	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(x1);
	fmid=(*func)(x2);
	if (f*fmid >= 0.0) NR_error("Root must be bracketed for bisection in NR_rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	NR_error("Too many bisections in NR_rtbis");
	return 0.0;
}
#undef JMAX
