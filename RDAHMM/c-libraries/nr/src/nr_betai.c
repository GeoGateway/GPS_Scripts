#include <math.h>

float NR_betai(float a, float b, float x)
{
	float NR_betacf(float a, float b, float x);
	float NR_gammln(float xx);
	void NR_error(char error_text[]);
	float bt;

	if (x < 0.0 || x > 1.0) NR_error("Bad x in routine NR_betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(NR_gammln(a+b)-NR_gammln(a)-NR_gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*NR_betacf(a,b,x)/a;
	else
		return 1.0-bt*NR_betacf(b,a,1.0-x)/b;
}
