#include <math.h>
#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float NR_qromo(float (*func)(float), float a, float b,
	float (*choose)(float(*)(float), float, float, int))
{
	void NR_polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void NR_error(char error_text[]);
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			NR_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	NR_error("Too many steps in routing NR_qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
