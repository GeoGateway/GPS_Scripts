#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr_util.h"
#define KBF 128

void NR_fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign)
{
	void NR_fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
	unsigned long j,j12,jk,k,kk,n=1,mm,kc=0,kd,ks,kr,nr,ns,nv;
	int cc,na,nb,nc,nd;
	float tempr,tempi,*afa,*afb,*afc;
	double wr,wi,wpr,wpi,wtemp,theta;
	static int mate[5] = {0,2,1,4,3};

	afa=NR_vector(1,KBF);
	afb=NR_vector(1,KBF);
	afc=NR_vector(1,KBF);
	for (j=1;j<=ndim;j++) {
		n *= nn[j];
		if (nn[j] <= 1) NR_error("invalid float or wrong ndim in NR_fourfs");
	}
	nv=1;
	jk=nn[nv];
	mm=n;
	ns=n/KBF;
	nr=ns >> 1;
	kd=KBF >> 1;
	ks=n;
	NR_fourew(file,&na,&nb,&nc,&nd);
	for (;;) {
		theta=isign*3.141592653589793/(n/mm);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		mm >>= 1;
		for (j12=1;j12<=2;j12++) {
			kr=0;
			do {
				cc=fread(&afa[1],sizeof(float),KBF,file[na]);
				if (cc != KBF) NR_error("read error in NR_fourfs");
				cc=fread(&afb[1],sizeof(float),KBF,file[nb]);
				if (cc != KBF) NR_error("read error in NR_fourfs");
				for (j=1;j<=KBF;j+=2) {
					tempr=((float)wr)*afb[j]-((float)wi)*afb[j+1];
					tempi=((float)wi)*afb[j]+((float)wr)*afb[j+1];
					afb[j]=afa[j]-tempr;
					afa[j] += tempr;
					afb[j+1]=afa[j+1]-tempi;
					afa[j+1] += tempi;
				}
				kc += kd;
				if (kc == mm) {
					kc=0;
					wr=(wtemp=wr)*wpr-wi*wpi+wr;
					wi=wi*wpr+wtemp*wpi+wi;
				}
				cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]);
				if (cc != KBF) NR_error("write error in NR_fourfs");
				cc=fwrite(&afb[1],sizeof(float),KBF,file[nd]);
				if (cc != KBF) NR_error("write error in NR_fourfs");
			} while (++kr < nr);
			if (j12 == 1 && ks != n && ks == KBF) {
				na=mate[na];
				nb=na;
			}
			if (nr == 0) break;
		}
		NR_fourew(file,&na,&nb,&nc,&nd);
		jk >>= 1;
		while (jk == 1) {
			mm=n;
			jk=nn[++nv];
		}
		ks >>= 1;
		if (ks > KBF) {
			for (j12=1;j12<=2;j12++) {
				for (kr=1;kr<=ns;kr+=ks/KBF) {
					for (k=1;k<=ks;k+=KBF) {
						cc=fread(&afa[1],sizeof(float),KBF,file[na]);
						if (cc != KBF) NR_error("read error in NR_fourfs");
						cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]);
						if (cc != KBF) NR_error("write error in NR_fourfs");
					}
					nc=mate[nc];
				}
				na=mate[na];
			}
			NR_fourew(file,&na,&nb,&nc,&nd);
		} else if (ks == KBF) nb=na;
		else break;
	}
	j=1;
	for (;;) {
		theta=isign*3.141592653589793/(n/mm);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		mm >>= 1;
		ks=kd;
		kd >>= 1;
		for (j12=1;j12<=2;j12++) {
			for (kr=1;kr<=ns;kr++) {
				cc=fread(&afc[1],sizeof(float),KBF,file[na]);
				if (cc != KBF) NR_error("read error in NR_fourfs");
				kk=1;
				k=ks+1;
				for (;;) {
					tempr=((float)wr)*afc[kk+ks]-((float)wi)*afc[kk+ks+1];
					tempi=((float)wi)*afc[kk+ks]+((float)wr)*afc[kk+ks+1];
					afa[j]=afc[kk]+tempr;
					afb[j]=afc[kk]-tempr;
					afa[++j]=afc[++kk]+tempi;
					afb[j++]=afc[kk++]-tempi;
					if (kk < k) continue;
					kc += kd;
					if (kc == mm) {
						kc=0;
						wr=(wtemp=wr)*wpr-wi*wpi+wr;
						wi=wi*wpr+wtemp*wpi+wi;
					}
					kk += ks;
					if (kk > KBF) break;
					else k=kk+ks;
				}
				if (j > KBF) {
					cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]);
					if (cc != KBF) NR_error("write error in NR_fourfs");
					cc=fwrite(&afb[1],sizeof(float),KBF,file[nd]);
					if (cc != KBF) NR_error("write error in NR_fourfs");
					j=1;
				}
			}
			na=mate[na];
		}
		NR_fourew(file,&na,&nb,&nc,&nd);
		jk >>= 1;
		if (jk > 1) continue;
		mm=n;
		do {
			if (nv < ndim) jk=nn[++nv];
			else {
				NR_free_vector(afc,1,KBF);
				NR_free_vector(afb,1,KBF);
				NR_free_vector(afa,1,KBF);
				return;
			}
		} while (jk == 1);
	}
}
#undef KBF
#undef NRANSI
