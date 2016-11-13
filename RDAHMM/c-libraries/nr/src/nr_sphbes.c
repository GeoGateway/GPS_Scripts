#include <math.h>
#define RTPIO2 1.2533141

void NR_sphbes(int n, float x, float *sj, float *sy, float *sjp, float *syp)
{
	void NR_bessjy(float x, float xnu, float *NR_rj, float *ry, float *NR_rjp,
		float *ryp);
	void NR_error(char error_text[]);
	float factor,order,NR_rj,NR_rjp,ry,ryp;

	if (n < 0 || x <= 0.0) NR_error("bad arguments in NR_sphbes");
	order=n+0.5;
	NR_bessjy(x,order,&NR_rj,&ry,&NR_rjp,&ryp);
	factor=RTPIO2/sqrt(x);
	*sj=factor*NR_rj;
	*sy=factor*ry;
	*sjp=factor*NR_rjp-(*sj)/(2.0*x);
	*syp=factor*ryp-(*sy)/(2.0*x);
}
#undef RTPIO2
