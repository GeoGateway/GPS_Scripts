#define NRANSI
#include "nr_util.h"

void NR_bcuint(float y[], float y1[], float y2[], float y12[], float x1l,
	float x1u, float x2l, float x2u, float x1, float x2, float *ansy,
	float *ansy1, float *ansy2)
{
	void NR_bcucof(float y[], float y1[], float y2[], float y12[], float d1,
		float d2, float **c);
	int i;
	float t,u,d1,d2,**c;

	c=NR_matrix(1,4,1,4);
	d1=x1u-x1l;
	d2=x2u-x2l;
	NR_bcucof(y,y1,y2,y12,d1,d2,c);
	if (x1u == x1l || x2u == x2l) NR_error("Bad input in routine NR_bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
		*ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
	}
	*ansy1 /= d1;
	*ansy2 /= d2;
	NR_free_matrix(c,1,4,1,4);
}
#undef NRANSI
