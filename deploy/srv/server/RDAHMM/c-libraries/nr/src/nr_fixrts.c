#include <math.h>
#include "nr_complex.h"
#define NMAX 100
#define ZERO NR_Complex(0.0,0.0)
#define ONE NR_Complex(1.0,0.0)

void NR_fixrts(float d[], int m)
{
	void NR_zroots(nr_fcomplex a[], int m, nr_fcomplex roots[], int polish);
	int i,j,polish;
	nr_fcomplex a[NMAX],roots[NMAX];

	a[m]=ONE;
	for (j=m-1;j>=0;j--)
		a[j]=NR_Complex(-d[m-j],0.0);
	polish=1;
	NR_zroots(a,m,roots,polish);
	for (j=1;j<=m;j++)
		if (NR_Cabs(roots[j]) > 1.0)
			roots[j]=NR_Cdiv(ONE,NR_Conjg(roots[j]));
	a[0]=NR_Csub(ZERO,roots[1]);
	a[1]=ONE;
	for (j=2;j<=m;j++) {
		a[j]=ONE;
		for (i=j;i>=2;i--)
			a[i-1]=NR_Csub(a[i-2],NR_Cmul(roots[j],a[i-1]));
		a[0]=NR_Csub(ZERO,NR_Cmul(roots[j],a[0]));
	}
	for (j=0;j<=m-1;j++)
		d[m-j] = -a[j].r;
}
#undef NMAX
#undef ZERO
#undef ONE
