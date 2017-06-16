#define NRANSI
#include "nr_util.h"
#define NPRE 1
#define NPOST 1
#define NGMAX 15

void NR_mglin(double **u, int n, int ncycle)
{
	void NR_addint(double **uf, double **uc, double **res, int nf);
	void NR_copy(double **aout, double **ain, int n);
	void NR_fill0(double **u, int n);
	void NR_interp(double **uf, double **uc, int nf);
	void NR_relax(double **u, double **rhs, int n);
	void NR_resid(double **res, double **u, double **rhs, int n);
	void NR_rstrct(double **uc, double **uf, int nc);
	void NR_slvsml(double **u, double **rhs);
	unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
	double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) NR_error("n-1 must be a power of 2 in NR_mglin.");
	if (ng > NGMAX) NR_error("increase NGMAX in NR_mglin.");
	nn=n/2+1;
	ngrid=ng-1;
	irho[ngrid]=NR_dmatrix(1,nn,1,nn);
	NR_rstrct(irho[ngrid],u,nn);
	while (nn > 3) {
		nn=nn/2+1;
		irho[--ngrid]=NR_dmatrix(1,nn,1,nn);
		NR_rstrct(irho[ngrid],irho[ngrid+1],nn);
	}
	nn=3;
	iu[1]=NR_dmatrix(1,nn,1,nn);
	irhs[1]=NR_dmatrix(1,nn,1,nn);
	NR_slvsml(iu[1],irho[1]);
	NR_free_dmatrix(irho[1],1,nn,1,nn);
	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=NR_dmatrix(1,nn,1,nn);
		irhs[j]=NR_dmatrix(1,nn,1,nn);
		ires[j]=NR_dmatrix(1,nn,1,nn);
		NR_interp(iu[j],iu[j-1],nn);
		NR_copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
		for (jcycle=1;jcycle<=ncycle;jcycle++) {
			nf=nn;
			for (jj=j;jj>=2;jj--) {
			for (jpre=1;jpre<=NPRE;jpre++)
				NR_relax(iu[jj],irhs[jj],nf);
			NR_resid(ires[jj],iu[jj],irhs[jj],nf);
			nf=nf/2+1;
			NR_rstrct(irhs[jj-1],ires[jj],nf);
			NR_fill0(iu[jj-1],nf);
			}
			NR_slvsml(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
			nf=2*nf-1;
			NR_addint(iu[jj],iu[jj-1],ires[jj],nf);
			for (jpost=1;jpost<=NPOST;jpost++)
				NR_relax(iu[jj],irhs[jj],nf);
			}
		}
	}
	NR_copy(u,iu[ngrid],n);
	for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) {
		NR_free_dmatrix(ires[j],1,nn,1,nn);
		NR_free_dmatrix(irhs[j],1,nn,1,nn);
		NR_free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng) NR_free_dmatrix(irho[j],1,nn,1,nn);
	}
	NR_free_dmatrix(irhs[1],1,3,1,3);
	NR_free_dmatrix(iu[1],1,3,1,3);
}
#undef NPRE
#undef NPOST
#undef NGMAX
#undef NRANSI
