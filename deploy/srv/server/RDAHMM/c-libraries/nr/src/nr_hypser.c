#include "nr_complex.h"
#define ONE NR_Complex(1.0,0.0)

void NR_hypser(nr_fcomplex a, nr_fcomplex b, nr_fcomplex c, nr_fcomplex z, nr_fcomplex *series,
	nr_fcomplex *deriv)
{
	void NR_error(char error_text[]);
	int n;
	nr_fcomplex aa,bb,cc,fac,temp;

	deriv->r=0.0;
	deriv->i=0.0;
	fac=NR_Complex(1.0,0.0);
	temp=fac;
	aa=a;
	bb=b;
	cc=c;
	for (n=1;n<=1000;n++) {
		fac=NR_Cmul(fac,NR_Cdiv(NR_Cmul(aa,bb),cc));
		deriv->r+=fac.r;
		deriv->i+=fac.i;
		fac=NR_Cmul(fac,NR_RCmul(1.0/n,z));
		*series=NR_Cadd(temp,fac);
		if (series->r == temp.r && series->i == temp.i) return;
		temp= *series;
		aa=NR_Cadd(aa,ONE);
		bb=NR_Cadd(bb,ONE);
		cc=NR_Cadd(cc,ONE);

	}
	NR_error("convergence failure in NR_hypser");
}
#undef ONE
