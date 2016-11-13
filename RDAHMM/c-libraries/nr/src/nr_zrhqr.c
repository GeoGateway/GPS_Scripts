#define NRANSI
#include "nr_util.h"
#define MAXM 50

void NR_zrhqr(float a[], int m, float rtr[], float rti[])
{
	void NR_balanc(float **a, int n);
	void NR_hqr(float **a, int n, float wr[], float wi[]);
	int j,k;
	float **hess,xr,xi;

	hess=NR_matrix(1,MAXM,1,MAXM);
	if (m > MAXM || a[m] == 0.0) NR_error("bad args in zrNR_hqr");
	for (k=1;k<=m;k++) {
		hess[1][k] = -a[m-k]/a[m];
		for (j=2;j<=m;j++) hess[j][k]=0.0;
		if (k != m) hess[k+1][k]=1.0;
	}
	NR_balanc(hess,m);
	NR_hqr(hess,m,rtr,rti);
	for (j=2;j<=m;j++) {
		xr=rtr[j];
		xi=rti[j];
		for (k=j-1;k>=1;k--) {
			if (rtr[k] <= xr) break;
			rtr[k+1]=rtr[k];
			rti[k+1]=rti[k];
		}
		rtr[k+1]=xr;
		rti[k+1]=xi;
	}
	NR_free_matrix(hess,1,MAXM,1,MAXM);
}
#undef MAXM
#undef NRANSI
