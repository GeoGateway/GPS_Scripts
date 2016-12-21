#include <math.h>
#define JMAX 20

float rtNR_newt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	void NR_error(char error_text[]);
	int j;
	float df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		(*funcd)(rtn,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			NR_error("Jumped out of brackets in rtNR_newt");
		if (fabs(dx) < xacc) return rtn;
	}
	NR_error("Maximum number of iterations exceeded in rtNR_newt");
	return 0.0;
}
#undef JMAX
