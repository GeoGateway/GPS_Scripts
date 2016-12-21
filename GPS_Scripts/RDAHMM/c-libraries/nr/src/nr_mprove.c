#define NRANSI
#include "nr_util.h"

void NR_mprove(float **a, float **alud, int n, int indx[], float b[], float x[])
{
	void NR_lubksb(float **a, int n, int *indx, float b[]);
	int j,i;
	double sdp;
	float *r;

	r=NR_vector(1,n);
	for (i=1;i<=n;i++) {
		sdp = -b[i];
		for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
		r[i]=sdp;
	}
	NR_lubksb(alud,n,indx,r);
	for (i=1;i<=n;i++) x[i] -= r[i];
	NR_free_vector(r,1,n);
}
#undef NRANSI
