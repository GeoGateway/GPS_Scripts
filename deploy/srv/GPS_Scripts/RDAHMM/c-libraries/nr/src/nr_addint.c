void NR_addint(double **uf, double **uc, double **res, int nf)
{
	void NR_interp(double **uf, double **uc, int nf);
	int i,j;

	NR_interp(res,uc,nf);
	for (j=1;j<=nf;j++)
		for (i=1;i<=nf;i++)
			uf[i][j] += res[i][j];
}
