#include <math.h>
#include "nr_complex.h"
#define EPS 2.0e-6
#define MAXM 100

void NR_zroots(nr_fcomplex a[], int m, nr_fcomplex roots[], int polish)
{
	void NR_laguer(nr_fcomplex a[], int m, nr_fcomplex *x, int *its);
	int i,its,j,jj;
	nr_fcomplex x,b,c,ad[MAXM];

	for (j=0;j<=m;j++) ad[j]=a[j];
	for (j=m;j>=1;j--) {
		x=NR_Complex(0.0,0.0);
		NR_laguer(ad,j,&x,&its);
		if (fabs(x.i) <= 2.0*EPS*fabs(x.r)) x.i=0.0;
		roots[j]=x;
		b=ad[j];
		for (jj=j-1;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=NR_Cadd(NR_Cmul(x,b),c);
		}
	}
	if (polish)
		for (j=1;j<=m;j++)
			NR_laguer(a,m,&roots[j],&its);
	for (j=2;j<=m;j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].r <= x.r) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}
#undef EPS
#undef MAXM
