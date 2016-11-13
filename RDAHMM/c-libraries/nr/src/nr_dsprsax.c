void NR_dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void NR_error(char error_text[]);
	unsigned long i,k;

	if (ija[1] != n+2) NR_error("NR_dsprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}
