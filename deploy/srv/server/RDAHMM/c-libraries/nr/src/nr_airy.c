#include <math.h>
#ifndef PI
#ifdef M_PI
#define PI  M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif
#define THIRD (1.0/3.0)
#define TWOTHR (2.0*THIRD)
#define ONOVRT 0.57735027

void NR_airy(float x, float *ai, float *bi, float *aip, float *bip)
{
	void NR_bessik(float x, float xnu, float *ri, float *rk, float *rip,
		float *rkp);
	void NR_bessjy(float x, float xnu, float *NR_rj, float *ry, float *NR_rjp,
		float *ryp);
	float absx,ri,rip,NR_rj,NR_rjp,rk,rkp,rootx,ry,ryp,z;

	absx=fabs(x);
	rootx=sqrt(absx);
	z=TWOTHR*absx*rootx;
	if (x > 0.0) {
		NR_bessik(z,THIRD,&ri,&rk,&rip,&rkp);
		*ai=rootx*ONOVRT*rk/PI;
		*bi=rootx*(rk/PI+2.0*ONOVRT*ri);
		NR_bessik(z,TWOTHR,&ri,&rk,&rip,&rkp);
		*aip = -x*ONOVRT*rk/PI;
		*bip=x*(rk/PI+2.0*ONOVRT*ri);
	} else if (x < 0.0) {
		NR_bessjy(z,THIRD,&NR_rj,&ry,&NR_rjp,&ryp);
		*ai=0.5*rootx*(NR_rj-ONOVRT*ry);
		*bi = -0.5*rootx*(ry+ONOVRT*NR_rj);
		NR_bessjy(z,TWOTHR,&NR_rj,&ry,&NR_rjp,&ryp);
		*aip=0.5*absx*(ONOVRT*ry+NR_rj);
		*bip=0.5*absx*(ONOVRT*NR_rj-ry);
	} else {
		*ai=0.35502805;
		*bi=(*ai)/ONOVRT;
		*aip = -0.25881940;
		*bip = -(*aip)/ONOVRT;
	}
}
#undef PI
#undef THIRD
#undef TWOTHR
#undef ONOVRT
