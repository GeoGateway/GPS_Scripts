#include <math.h>
#define NRANSI
#include "nr_util.h"
#define TINY 1.0e-30

void NR_cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy,
	float *hygx, float *hxgy, float *uygx, float *uxgy, float *uxy)
{
	int i,j;
	float sum=0.0,p,*sumi,*sumj;

	sumi=NR_vector(1,ni);
	sumj=NR_vector(1,nj);
	for (i=1;i<=ni;i++) {
		sumi[i]=0.0;
		for (j=1;j<=nj;j++) {
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
	}
	for (j=1;j<=nj;j++) {
		sumj[j]=0.0;
		for (i=1;i<=ni;i++)
			sumj[j] += nn[i][j];
	}
	*hx=0.0;
	for (i=1;i<=ni;i++)
		if (sumi[i]) {
			p=sumi[i]/sum;
			*hx -= p*log(p);
		}
	*hy=0.0;
	for (j=1;j<=nj;j++)
		if (sumj[j]) {
			p=sumj[j]/sum;
			*hy -= p*log(p);
		}
	*h=0.0;
	for (i=1;i<=ni;i++)
		for (j=1;j<=nj;j++)
			if (nn[i][j]) {
				p=nn[i][j]/sum;
				*h -= p*log(p);
			}
	*hygx=(*h)-(*hx);
	*hxgy=(*h)-(*hy);
	*uygx=(*hy-*hygx)/(*hy+TINY);
	*uxgy=(*hx-*hxgy)/(*hx+TINY);
	*uxy=2.0*(*hx+*hy-*h)/(*hx+*hy+TINY);
	NR_free_vector(sumj,1,nj);
	NR_free_vector(sumi,1,ni);
}
#undef TINY
#undef NRANSI
