#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TFACTR 0.9
#define ALEN(a,b,c,d) sqrt(((b)-(a))*((b)-(a))+((d)-(c))*((d)-(c)))

void NR_anneal(float x[], float y[], int iorder[], int ncity)
{
	int NR_irbit1(unsigned long *iseed);
	int NR_metrop(float de, float t);
	float NR_ran3(long *idum);
	float NR_revcst(float x[], float y[], int iorder[], int ncity, int n[]);
	void NR_reverse(int iorder[], int ncity, int n[]);
	float NR_trncst(float x[], float y[], int iorder[], int ncity, int n[]);
	void NR_trnspt(int iorder[], int ncity, int n[]);
	int ans,nover,nlimit,i1,i2;
	int i,j,k,nsucc,nn,idec;
	static int n[7];
	long idum;
	unsigned long iseed;
	float path,de,t;

	nover=100*ncity;
	nlimit=10*ncity;
	path=0.0;
	t=0.5;
	for (i=1;i<ncity;i++) {
		i1=iorder[i];
		i2=iorder[i+1];
		path += ALEN(x[i1],x[i2],y[i1],y[i2]);
	}
	i1=iorder[ncity];
	i2=iorder[1];
	path += ALEN(x[i1],x[i2],y[i1],y[i2]);
	idum = -1;
	iseed=111;
	for (j=1;j<=100;j++) {
		nsucc=0;
		for (k=1;k<=nover;k++) {
			do {
				n[1]=1+(int) (ncity*NR_ran3(&idum));
				n[2]=1+(int) ((ncity-1)*NR_ran3(&idum));
				if (n[2] >= n[1]) ++n[2];
				nn=1+((n[1]-n[2]+ncity-1) % ncity);
			} while (nn<3);
			idec=NR_irbit1(&iseed);
			if (idec == 0) {
				n[3]=n[2]+(int) (abs(nn-2)*NR_ran3(&idum))+1;
				n[3]=1+((n[3]-1) % ncity);
				de=NR_trncst(x,y,iorder,ncity,n);
				ans=NR_metrop(de,t);
				if (ans) {
					++nsucc;
					path += de;
					NR_trnspt(iorder,ncity,n);
				}
			} else {
				de=NR_revcst(x,y,iorder,ncity,n);
				ans=NR_metrop(de,t);
				if (ans) {
					++nsucc;
					path += de;
					NR_reverse(iorder,ncity,n);
				}
			}
			if (nsucc >= nlimit) break;
		}
		printf("\n %s %10.6f %s %12.6f \n","T =",t,
			"	 Path Length =",path);
		printf("Successful Moves: %6d\n",nsucc);
		t *= TFACTR;
		if (nsucc == 0) return;
	}
}
#undef TFACTR
#undef ALEN
