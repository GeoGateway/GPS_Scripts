#include <math.h>
#define NRANSI
#include "nr_util.h"
#define BIG 1.0e30

void NR_pade(double cof[], int n, float *NR_resid)
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	void NR_ludcmp(float **a, int n, int *indx, float *d);
	void NR_mprove(float **a, float **alud, int n, int indx[], float b[],
		float x[]);
	int j,k,*indx;
	float d,rr,rrold,sum,**q,**qlu,*x,*y,*z;

	indx=NR_ivector(1,n);
	q=NR_matrix(1,n,1,n);
	qlu=NR_matrix(1,n,1,n);
	x=NR_vector(1,n);
	y=NR_vector(1,n);
	z=NR_vector(1,n);
	for (j=1;j<=n;j++) {
		y[j]=x[j]=cof[n+j];
		for (k=1;k<=n;k++) {
			q[j][k]=cof[j-k+n];
			qlu[j][k]=q[j][k];
		}
	}
	NR_ludcmp(qlu,n,indx,&d);
	NR_lubksb(qlu,n,indx,x);
	rr=BIG;
	do {
		rrold=rr;
		for (j=1;j<=n;j++) z[j]=x[j];
		NR_mprove(q,qlu,n,indx,y,x);
		for (rr=0.0,j=1;j<=n;j++)
			rr += NR_SQR(z[j]-x[j]);
	} while (rr < rrold);
	*NR_resid=sqrt(rr);
	for (k=1;k<=n;k++) {
		for (sum=cof[k],j=1;j<=k;j++) sum -= x[j]*cof[k-j];
		y[k]=sum;
	}
	for (j=1;j<=n;j++) {
		cof[j]=y[j];
		cof[j+n] = -x[j];
	}
	NR_free_vector(z,1,n);
	NR_free_vector(y,1,n);
	NR_free_vector(x,1,n);
	NR_free_matrix(qlu,1,n,1,n);
	NR_free_matrix(q,1,n,1,n);
	NR_free_ivector(indx,1,n);
}
#undef BIG
#undef NRANSI
