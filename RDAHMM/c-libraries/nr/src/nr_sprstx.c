void NR_sprstx(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n)
{
	void NR_error(char error_text[]);
	unsigned long i,j,k;

	if (ija[1] != n+2) NR_error("mismatched vector and matrix in NR_sprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}
