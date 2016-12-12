/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file complex.c.  Do not confuse this file with the same-named
   file complex.c that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to NR_select the
   correct version.  *This* file contains only ANSI C.               */

#include <math.h>

typedef struct NR_FCOMPLEX {float r,i;} nr_fcomplex;

nr_fcomplex NR_Cadd(nr_fcomplex a, nr_fcomplex b)
{
	nr_fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

nr_fcomplex NR_Csub(nr_fcomplex a, nr_fcomplex b)
{
	nr_fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


nr_fcomplex NR_Cmul(nr_fcomplex a, nr_fcomplex b)
{
	nr_fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

nr_fcomplex NR_Complex(float re, float im)
{
	nr_fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

nr_fcomplex NR_Conjg(nr_fcomplex z)
{
	nr_fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

nr_fcomplex NR_Cdiv(nr_fcomplex a, nr_fcomplex b)
{
	nr_fcomplex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

float NR_Cabs(nr_fcomplex z)
{
	float x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

nr_fcomplex NR_Csqrt(nr_fcomplex z)
{
	nr_fcomplex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

nr_fcomplex NR_RCmul(float x, nr_fcomplex a)
{
	nr_fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}
