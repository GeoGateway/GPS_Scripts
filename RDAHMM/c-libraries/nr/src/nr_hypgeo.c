#include <math.h>
#include "nr_complex.h"
#define NRANSI
#include "nr_util.h"
#define EPS 1.0e-6

nr_fcomplex aa,bb,cc,z0,dz;

int kmax,kount;
float *xp,**yp,dxsav;

nr_fcomplex NR_hypgeo(nr_fcomplex a, nr_fcomplex b, nr_fcomplex c, nr_fcomplex z)
{
	void NR_bsstep(float y[], float dydx[], int nv, float *xx, float htry,
		float eps, float yscal[], float *hdid, float *hnext,
		void (*derivs)(float, float [], float []));
	void NR_hypdrv(float s, float yy[], float dyyds[]);
	void NR_hypser(nr_fcomplex a, nr_fcomplex b, nr_fcomplex c, nr_fcomplex z,
		nr_fcomplex *series, nr_fcomplex *deriv);
	void NR_odeint(float ystart[], int nvar, float x1, float x2,
		float eps, float h1, float hmin, int *nok, int *nbad,
		void (*derivs)(float, float [], float []),
		void (*NR_rkqs)(float [], float [], int, float *, float, float,
		float [], float *, float *, void (*)(float, float [], float [])));
	int nbad,nok;
	nr_fcomplex ans,y[3];
	float *yy;

	kmax=0;
	if (z.r*z.r+z.i*z.i <= 0.25) {
		NR_hypser(a,b,c,z,&ans,&y[2]);
		return ans;
	}
	else if (z.r < 0.0) z0=NR_Complex(-0.5,0.0);
	else if (z.r <= 1.0) z0=NR_Complex(0.5,0.0);
	else z0=NR_Complex(0.0,z.i >= 0.0 ? 0.5 : -0.5);
	aa=a;
	bb=b;
	cc=c;
	dz=NR_Csub(z,z0);
	NR_hypser(aa,bb,cc,z0,&y[1],&y[2]);
	yy=NR_vector(1,4);
	yy[1]=y[1].r;
	yy[2]=y[1].i;
	yy[3]=y[2].r;
	yy[4]=y[2].i;
	NR_odeint(yy,4,0.0,1.0,EPS,0.1,0.0001,&nok,&nbad,NR_hypdrv,NR_bsstep);
	y[1]=NR_Complex(yy[1],yy[2]);
	NR_free_vector(yy,1,4);
	return y[1];
}
#undef EPS
#undef NRANSI
