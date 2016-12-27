extern unsigned long ija[];
extern double sa[];

void NR_atimes(unsigned long n, double x[], double r[], int itrnsp)
{
	void NR_dsprsax(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);
	void NR_dsprstx(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);

	if (itrnsp) NR_dsprstx(sa,ija,x,r,n);
	else NR_dsprsax(sa,ija,x,r,n);
}
