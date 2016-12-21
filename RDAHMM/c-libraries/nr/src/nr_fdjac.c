#include <math.h>
#define NRANSI
#include "nr_util.h"
#define EPS 1.0e-4

void NR_fdjac(int n, float x[], float fvec[], float **df,
	void (*vecfunc)(int, float [], float []))
{
	int i,j;
	float h,temp,*f;

	f=NR_vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h=EPS;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	NR_free_vector(f,1,n);
}
#undef EPS
#undef NRANSI
