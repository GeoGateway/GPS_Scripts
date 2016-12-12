#include <stdlib.h>
#include <math.h>
#define NRANSI
#include "nr_util.h"
#define PFAC 0.1
#define MNPT 15
#define MNBS 60
#define TINY 1.0e-30
#define BIG 1.0e30

static long iran=0;

void NR_miser(float (*func)(float []), float regn[], int ndim, unsigned long npts,
	float dith, float *ave, float *var)
{
	void NR_ranpt(float pt[], float regn[], int n);
	float *regn_temp;
	unsigned long n,npre,nptl,nptr;
	int j,jb;
	float avel,varl;
	float fracl,fval;
	float rgl,rgm,rgr,s,sigl,siglb,sigr,sigrb;
	float sum,sumb,summ,summ2;
	float *fmaxl,*fmaxr,*NR_fminl,*NR_fminr;
	float *pt,*rmid;

	pt=NR_vector(1,ndim);
	if (npts < MNBS) {
		summ=summ2=0.0;
		for (n=1;n<=npts;n++) {
			NR_ranpt(pt,regn,ndim);
			fval=(*func)(pt);
			summ += fval;
			summ2 += fval * fval;
		}
		*ave=summ/npts;
		*var=NR_FMAX(TINY,(summ2-summ*summ/npts)/(npts*npts));
	}
	else {
		rmid=NR_vector(1,ndim);
		npre=NR_LMAX((unsigned long)(npts*PFAC),MNPT);
		fmaxl=NR_vector(1,ndim);
		fmaxr=NR_vector(1,ndim);
		NR_fminl=NR_vector(1,ndim);
		NR_fminr=NR_vector(1,ndim);
		for (j=1;j<=ndim;j++) {
			iran=(iran*2661+36979) % 175000;
			s=NR_SIGN(dith,(float)(iran-87500));
			rmid[j]=(0.5+s)*regn[j]+(0.5-s)*regn[ndim+j];
			NR_fminl[j]=NR_fminr[j]=BIG;
			fmaxl[j]=fmaxr[j] = -BIG;
		}
		for (n=1;n<=npre;n++) {
			NR_ranpt(pt,regn,ndim);
			fval=(*func)(pt);
			for (j=1;j<=ndim;j++) {
				if (pt[j]<=rmid[j]) {
					NR_fminl[j]=NR_FMIN(NR_fminl[j],fval);
					fmaxl[j]=NR_FMAX(fmaxl[j],fval);
				}
				else {
					NR_fminr[j]=NR_FMIN(NR_fminr[j],fval);
					fmaxr[j]=NR_FMAX(fmaxr[j],fval);
				}
			}
		}
		sumb=BIG;
		jb=0;
		siglb=sigrb=1.0;
		for (j=1;j<=ndim;j++) {
			if (fmaxl[j] > NR_fminl[j] && fmaxr[j] > NR_fminr[j]) {
				sigl=NR_FMAX(TINY,pow(fmaxl[j]-NR_fminl[j],2.0/3.0));
				sigr=NR_FMAX(TINY,pow(fmaxr[j]-NR_fminr[j],2.0/3.0));
				sum=sigl+sigr;
				if (sum<=sumb) {
					sumb=sum;
					jb=j;
					siglb=sigl;
					sigrb=sigr;
				}
			}
		}
		NR_free_vector(NR_fminr,1,ndim);
		NR_free_vector(NR_fminl,1,ndim);
		NR_free_vector(fmaxr,1,ndim);
		NR_free_vector(fmaxl,1,ndim);
		if (!jb) jb=1+(ndim*iran)/175000;
		rgl=regn[jb];
		rgm=rmid[jb];
		rgr=regn[ndim+jb];
		fracl=fabs((rgm-rgl)/(rgr-rgl));
		nptl=(unsigned long)(MNPT+(npts-npre-2*MNPT)*fracl*siglb
			/(fracl*siglb+(1.0-fracl)*sigrb));
		nptr=npts-npre-nptl;
		regn_temp=NR_vector(1,2*ndim);
		for (j=1;j<=ndim;j++) {
			regn_temp[j]=regn[j];
			regn_temp[ndim+j]=regn[ndim+j];
		}
		regn_temp[ndim+jb]=rmid[jb];
		NR_miser(func,regn_temp,ndim,nptl,dith,&avel,&varl);
		regn_temp[jb]=rmid[jb];
		regn_temp[ndim+jb]=regn[ndim+jb];
		NR_miser(func,regn_temp,ndim,nptr,dith,ave,var);
		NR_free_vector(regn_temp,1,2*ndim);
		*ave=fracl*avel+(1-fracl)*(*ave);
		*var=fracl*fracl*varl+(1-fracl)*(1-fracl)*(*var);
		NR_free_vector(rmid,1,ndim);
	}
	NR_free_vector(pt,1,ndim);
}
#undef MNPT
#undef MNBS
#undef TINY
#undef BIG
#undef PFAC
#undef NRANSI
