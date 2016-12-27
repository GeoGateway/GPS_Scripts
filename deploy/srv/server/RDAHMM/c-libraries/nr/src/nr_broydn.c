#include <math.h>
#define NRANSI
#include "nr_util.h"
#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TONR_LMIN 1.0e-6
#define FREERETURN {NR_free_vector(fvec,1,n);NR_free_vector(xold,1,n);\
	NR_free_vector(w,1,n);NR_free_vector(t,1,n);NR_free_vector(s,1,n);\
	NR_free_matrix(r,1,n,1,n);NR_free_matrix(qt,1,n,1,n);NR_free_vector(p,1,n);\
	NR_free_vector(g,1,n);NR_free_vector(fvcold,1,n);NR_free_vector(d,1,n);\
	NR_free_vector(c,1,n);return;}

int nn;
float *fvec;
void (*nrfuncv)(int n, float v[], float f[]);

void NR_broydn(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []))
{
	void NR_fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float NR_fmin(float x[]);
	void NR_lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void NR_qrdcmp(float **a, int n, float *c, float *d, int *sing);
	void NR_qrupdt(float **r, float **qt, int n, float u[], float v[]);
	void NR_rsolv(float **a, int n, float d[], float b[]);
	int i,its,j,k,restrt,sing,skip;
	float den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
	float *g,*p,**qt,**r,*s,*t,*w,*xold;

	c=NR_vector(1,n);
	d=NR_vector(1,n);
	fvcold=NR_vector(1,n);
	g=NR_vector(1,n);
	p=NR_vector(1,n);
	qt=NR_matrix(1,n,1,n);
	r=NR_matrix(1,n,1,n);
	s=NR_vector(1,n);
	t=NR_vector(1,n);
	w=NR_vector(1,n);
	xold=NR_vector(1,n);
	fvec=NR_vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=NR_fmin(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += NR_SQR(x[i]);
	stpmax=STPMX*NR_FMAX(sqrt(sum),(float)n);
	restrt=1;
	for (its=1;its<=MAXITS;its++) {
		if (restrt) {
			NR_fdjac(n,x,fvec,r,vecfunc);
			NR_qrdcmp(r,n,c,d,&sing);
			if (sing) NR_error("singular Jacobian in NR_broydn");
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) qt[i][j]=0.0;
				qt[i][i]=1.0;
			}
			for (k=1;k<n;k++) {
				if (c[k]) {
					for (j=1;j<=n;j++) {
						sum=0.0;
						for (i=k;i<=n;i++)
							sum += r[i][k]*qt[i][j];
						sum /= c[k];
						for (i=k;i<=n;i++)
							qt[i][j] -= sum*r[i][k];
					}
				}
			}
			for (i=1;i<=n;i++) {
				r[i][i]=d[i];
				for (j=1;j<i;j++) r[i][j]=0.0;
			}
		} else {
			for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
			for (i=1;i<=n;i++) {
				for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
				t[i]=sum;
			}
			skip=1;
			for (i=1;i<=n;i++) {
				for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
				w[i]=fvec[i]-fvcold[i]-sum;
				if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;
				else w[i]=0.0;
			}
			if (!skip) {
				for (i=1;i<=n;i++) {
					for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
					t[i]=sum;
				}
				for (den=0.0,i=1;i<=n;i++) den += NR_SQR(s[i]);
				for (i=1;i<=n;i++) s[i] /= den;
				NR_qrupdt(r,qt,n,t,s);
				for (i=1;i<=n;i++) {
					if (r[i][i] == 0.0) NR_error("r singular in NR_broydn");
					d[i]=r[i][i];
				}
			}
		}
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
			g[i]=sum;
		}
		for (i=n;i>=1;i--) {
			for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) {
			xold[i]=x[i];
			fvcold[i]=fvec[i];
		}
		fold=f;
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
			p[i] = -sum;
		}
		NR_rsolv(r,n,d,p);
		NR_lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,NR_fmin);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			if (restrt) FREERETURN
			else {
				test=0.0;
				den=NR_FMAX(f,0.5*n);
				for (i=1;i<=n;i++) {
					temp=fabs(g[i])*NR_FMAX(fabs(x[i]),1.0)/den;
					if (temp > test) test=temp;
				}
				if (test < TONR_LMIN) FREERETURN
				else restrt=1;
			}
		} else {
			restrt=0;
			test=0.0;
			for (i=1;i<=n;i++) {
				temp=(fabs(x[i]-xold[i]))/NR_FMAX(fabs(x[i]),1.0);
				if (temp > test) test=temp;
			}
			if (test < TOLX) FREERETURN
		}
	}
	NR_error("MAXITS exceeded in NR_broydn");
	FREERETURN
}
#undef MAXITS
#undef EPS
#undef TOLF
#undef TONR_LMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
