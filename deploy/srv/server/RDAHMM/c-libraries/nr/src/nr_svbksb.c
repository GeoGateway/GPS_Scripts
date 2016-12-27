#define NRANSI
#include "nr_util.h"

void NR_svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[])
{
	int jj,j,i;
	float s,*tmp;

	tmp=NR_vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	NR_free_vector(tmp,1,n);
}
#undef NRANSI
