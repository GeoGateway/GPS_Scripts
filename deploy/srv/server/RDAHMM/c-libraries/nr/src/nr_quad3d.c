static float xsav,ysav;
static float (*nrfunc)(float,float,float);

float NR_quad3d(float (*func)(float, float, float), float x1, float x2)
{
	float NR_qgaus(float (*func)(float), float a, float b);
	float f1(float x);

	nrfunc=func;
	return NR_qgaus(f1,x1,x2);
}

float f1(float x)
{
	float NR_qgaus(float (*func)(float), float a, float b);
	float f2(float y);
	float yy1(float),yy2(float);

	xsav=x;
	return NR_qgaus(f2,yy1(x),yy2(x));
}

float f2(float y)
{
	float NR_qgaus(float (*func)(float), float a, float b);
	float f3(float z);
	float z1(float,float),z2(float,float);

	ysav=y;
	return NR_qgaus(f3,z1(xsav,y),z2(xsav,y));
}

float f3(float z)
{
	return (*nrfunc)(xsav,ysav,z);
}
