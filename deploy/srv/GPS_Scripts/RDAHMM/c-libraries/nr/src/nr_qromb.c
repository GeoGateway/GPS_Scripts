#include <math.h>
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

float NR_qromb(float (*func)(float), float a, float b)
{
	void NR_polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	float NR_trapzd(float (*func)(float), float a, float b, int n);
	void NR_error(char error_text[]);
	float ss,dss;
	float s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=NR_trapzd(func,a,b,j);
		if (j >= K) {
			NR_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	NR_error("Too many steps in routine NR_qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
