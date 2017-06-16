#define NRANSI
#include "nr_util.h"
#define NPRE 1
#define NPOST 1
#define ALPHA 0.33
#define NGMAX 15

void NR_mgfas(double **u, int n, int maxcyc)
{
	double NR_anorm2(double **a, int n);
	void NR_copy(double **aout, double **ain, int n);
	void NR_interp(double **uf, double **uc, int nf);
	void NR_lop(double **out, double **u, int n);
	void NR_matadd(double **a, double **b, double **c, int n);
	void NR_matsub(double **a, double **b, double **c, int n);
	void NR_relax2(double **u, double **rhs, int n);
	void NR_rstrct(double **uc, double **uf, int nc);
	void NR_slvsm2(double **u, double **rhs);
	unsigned int j,jcycle,jj,jm1,jpost,jpre,nf,ng=0,ngrid,nn;
	double **irho[NGMAX+1],**irhs[NGMAX+1],**itau[NGMAX+1],
		**itemp[NGMAX+1],**iu[NGMAX+1];
	double res,trerr;

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) NR_error("n-1 must be a power of 2 in NR_mgfas.");
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
	itau[1]=NR_dmatrix(1,nn,1,nn);
	itemp[1]=NR_dmatrix(1,nn,1,nn);
	NR_slvsm2(iu[1],irho[1]);
	NR_free_dmatrix(irho[1],1,nn,1,nn);
	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=NR_dmatrix(1,nn,1,nn);
		irhs[j]=NR_dmatrix(1,nn,1,nn);
		itau[j]=NR_dmatrix(1,nn,1,nn);
		itemp[j]=NR_dmatrix(1,nn,1,nn);
		NR_interp(iu[j],iu[j-1],nn);
		NR_copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
		for (jcycle=1;jcycle<=maxcyc;jcycle++) {
		nf=nn;
			for (jj=j;jj>=2;jj--) {
				for (jpre=1;jpre<=NPRE;jpre++)
					NR_relax2(iu[jj],irhs[jj],nf);
				NR_lop(itemp[jj],iu[jj],nf);
				nf=nf/2+1;
				jm1=jj-1;
				NR_rstrct(itemp[jm1],itemp[jj],nf);
				NR_rstrct(iu[jm1],iu[jj],nf);
				NR_lop(itau[jm1],iu[jm1],nf);
				NR_matsub(itau[jm1],itemp[jm1],itau[jm1],nf);
				if (jj == j)
					trerr=ALPHA*NR_anorm2(itau[jm1],nf);
				NR_rstrct(irhs[jm1],irhs[jj],nf);
				NR_matadd(irhs[jm1],itau[jm1],irhs[jm1],nf);
			}
			NR_slvsm2(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
			jm1=jj-1;
			NR_rstrct(itemp[jm1],iu[jj],nf);
			NR_matsub(iu[jm1],itemp[jm1],itemp[jm1],nf);
			nf=2*nf-1;
			NR_interp(itau[jj],itemp[jm1],nf);
			NR_matadd(iu[jj],itau[jj],iu[jj],nf);
			for (jpost=1;jpost<=NPOST;jpost++)
				NR_relax2(iu[jj],irhs[jj],nf);
			}
			NR_lop(itemp[j],iu[j],nf);
			NR_matsub(itemp[j],irhs[j],itemp[j],nf);
			res=NR_anorm2(itemp[j],nf);
			if (res < trerr) break;
		}
	}
	NR_copy(u,iu[ngrid],n);
	for (nn=n,j=ng;j>=1;j--,nn=nn/2+1) {
		NR_free_dmatrix(itemp[j],1,nn,1,nn);
		NR_free_dmatrix(itau[j],1,nn,1,nn);
		NR_free_dmatrix(irhs[j],1,nn,1,nn);
		NR_free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng && j != 1) NR_free_dmatrix(irho[j],1,nn,1,nn);
	}
}
#undef NGMAX
#undef NPRE
#undef NPOST
#undef ALPHA
#undef NRANSI
