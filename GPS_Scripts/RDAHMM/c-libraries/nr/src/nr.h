/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nr.h.  Do not confuse this file with the same-named
   file nr.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_H_
#define _NR_H_

#ifndef lint
static char nr_hdr_rcsid[] = "$Id: nr.h,v 1.5 1997/07/29 15:25:49 granat Exp granat $";
#endif
/* 
 * $Log: nr.h,v $
 * Revision 1.5  1997/07/29 15:25:49  granat
 * edited to match new "NR_" paradigm
 *
 * Revision 1.4  1997/05/19 18:30:31  granat
 * changed da_select prototype to "select"
 *
 * Revision 1.3  1996/09/24 17:31:59  granat
 * Changed prototypes of hacked functions to have da_ prefix
 *
 * Revision 1.2  1996/07/25 18:14:42  granat
 * added prototypes for svdcmp2, svdcmp3
 *
 * Revision 1.1  1996/07/25 17:27:35  agray
 * Initial revision
 *
 * */

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>
#include "nr_complex.h"
#include "nr_util.h"

void NR_addint(double **uf, double **uc, double **res, int nf);
void NR_airy(float x, float *ai, float *bi, float *aip, float *bip);
void NR_amebsa(float **p, float y[], int ndim, float pb[],	float *yb,
	float ftol, float (*funk)(float []), int *iter, float temptr);
void NR_amoeba(float **p, float y[], int ndim, float ftol,
	float (*funk)(float []), int *iter);
float NR_amotry(float **p, float y[], float psum[], int ndim,
	float (*funk)(float []), int ihi, float fac);
float NR_amotsa(float **p, float y[], float psum[], int ndim, float pb[],
	float *yb, float (*funk)(float []), int ihi, float *yhi, float fac);
void NR_anneal(float x[], float y[], int iorder[], int ncity);
double NR_anorm2(double **a, int n);
void NR_arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
	arithcode *acode);
void NR_arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void NR_arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void NR_asolve(unsigned long n, double b[], double x[], int itrnsp);
void NR_atimes(unsigned long n, double x[], double r[], int itrnsp);
void NR_avevar(float data[], unsigned long n, float *ave, float *var);
void NR_balanc(float **a, int n);
void NR_banbks(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float b[]);
void NR_bandec(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float *d);
void NR_banmul(float **a, unsigned long n, int m1, int m2, float x[], float b[]);
void NR_bcucof(float y[], float y1[], float y2[], float y12[], float d1,
	float d2, float **c);
void NR_bcuint(float y[], float y1[], float y2[], float y12[],
	float x1l, float x1u, float x2l, float x2u, float x1,
	float x2, float *ansy, float *ansy1, float *ansy2);
void NR_beschb(double x, double *gam1, double *gam2, double *gampl,
	double *gammi);
float NR_bessi(int n, float x);
float NR_bessi0(float x);
float NR_bessi1(float x);
void NR_bessik(float x, float xnu, float *ri, float *rk, float *rip,
	float *rkp);
float NR_bessj(int n, float x);
float NR_bessj0(float x);
float NR_bessj1(float x);
void NR_bessjy(float x, float xnu, float *rj, float *ry, float *rjp,
	float *ryp);
float NR_bessk(int n, float x);
float NR_bessk0(float x);
float NR_bessk1(float x);
float NR_bessy(int n, float x);
float NR_bessy0(float x);
float NR_bessy1(float x);
float NR_beta(float z, float w);
float NR_betacf(float a, float b, float x);
float NR_betai(float a, float b, float x);
float NR_bico(int n, int k);
void NR_bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
float NR_bnldev(float pp, int n, long *idum);
float NR_brent(float ax, float bx, float cx,
	float (*f)(float), float tol, float *xmin);
void NR_broydn(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void NR_bsstep(float y[], float dydx[], int nv, float *xx, float htry,
	float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void NR_caldat(long julian, int *mm, int *id, int *iyyy);
void NR_chder(float a, float b, float c[], float cder[], int n);
float NR_chebev(float a, float b, float c[], int m, float x);
void NR_chebft(float a, float b, float c[], int n, float (*func)(float));
void NR_chebpc(float c[], float d[], int n);
void NR_chint(float a, float b, float c[], float cint[], int n);
float NR_chixy(float bang);
void NR_choldc(float **a, int n, float p[]);
void NR_cholsl(float **a, int n, float p[], float b[], float x[]);
void NR_chsone(float bins[], float ebins[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void NR_chstwo(float bins1[], float bins2[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void NR_cisi(float x, float *ci, float *si);
void NR_cntab1(int **nn, int ni, int nj, float *chisq,
	float *df, float *prob, float *cramrv, float *ccc);
void NR_cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy,
	float *hygx, float *hxgy, float *uygx, float *uxgy, float *uxy);
void NR_convlv(float data[], unsigned long n, float respns[], unsigned long m,
	int isign, float ans[]);
void NR_copy(double **aout, double **ain, int n);
void NR_correl(float data1[], float data2[], unsigned long n, float ans[]);
void NR_cosft(float y[], int n, int isign);
void NR_cosft1(float y[], int n);
void NR_cosft2(float y[], int n, int isign);
void NR_covsrt(float **covar, int ma, int ia[], int mfit);
void NR_crank(unsigned long n, float w[], float *s);
void NR_cyclic(float a[], float b[], float c[], float alpha, float beta,
	float r[], float x[], unsigned long n);
void NR_daub4(float a[], unsigned long n, int isign);
float NR_dawson(float x);
float NR_dbrent(float ax, float bx, float cx,
	float (*f)(float), float (*df)(float), float tol, float *xmin);
void NR_dcholdc(double **a, int n, double p[]);
void NR_dcholsl(double **a, int n, double p[], double b[], double x[]);
void NR_ddpoly(float c[], int nc, float x, float pd[], int nd);
int NR_decchk(char string[], int n, char *ch);
void NR_derivs(float x, float y[], float dydx[]);
float NR_df1dim(float x);
void NR_dfour1(double data[], unsigned long nn, int isign);
void NR_dfpmin(float p[], int n, float gtol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
float NR_dfridr(float (*func)(float), float x, float h, float *err);
void NR_dftcor(float w, float delta, float a, float b, float endpts[],
	float *corre, float *corim, float *corfac);
void NR_dftint(float (*func)(float), float a, float b, float w,
	float *cosint, float *sinint);
double NR_dgammln(double xx);
double NR_dgammp(double a, double x);
double NR_dgasdev(long *idum);
void NR_dgcf(double *gammcf, double a, double x, double *gln);
void NR_dgser(double *gamser, double a, double x, double *gln);
void NR_difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, float **s, float **y);
void NR_djacobi(double **a, int n, double d[], double **v, int *nrot);
void NR_dlinmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float[]));
void NR_dlubksb(double **a, int n, int *indx, double b[]);
void NR_dludcmp(double **a, int n, int *indx, double *d);
double NR_dpythag(double a, double b);
void NR_dqrdcmp(double **a, int n, double *c, double *d, int *sing);
void NR_drealft(double data[], unsigned long n, int isign);
void NR_dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void NR_dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void NR_dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void NR_dsvdcmp(double **a, int m, int n, double w[], double **v);
void NR_eclass(int nf[], int n, int lista[], int listb[], int m);
void NR_eclazz(int nf[], int n, int (*equiv)(int, int));
float NR_ei(float x);
void NR_eigsrt(float d[], float **v, int n);
float NR_elle(float phi, float ak);
float NR_ellf(float phi, float ak);
float NR_ellpi(float phi, float en, float ak);
void NR_elmhes(float **a, int n);
float NR_erfcc(float x);
float NR_erff(float x);
float NR_erffc(float x);
void NR_eulsum(float *sum, float term, int jterm, float wksp[]);
float NR_evlmem(float fdt, float d[], int m, float xms);
float NR_expdev(long *idum);
float NR_expint(int n, float x);
float NR_f1(float x);
float NR_f1dim(float x);
float NR_f2(float y);
float NR_f3(float z);
float NR_factln(int n);
float NR_factrl(int n);
void NR_fasper(float x[], float y[], unsigned long n, float ofac, float hifac,
	float wk1[], float wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, float *prob);
void NR_fdjac(int n, float x[], float fvec[], float **df,
	void (*vecfunc)(int, float [], float []));
void NR_fgauss(float x, float a[], float *y, float dyda[], int na);
void NR_fill0(double **u, int n);
void NR_fit(float x[], float y[], int ndata, float sig[], int mwt,
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void NR_fitexy(float x[], float y[], int ndat, float sigx[], float sigy[],
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void NR_fixrts(float d[], int m);
void NR_fleg(float x, float pl[], int nl);
void NR_flmoon(int n, int nph, long *jd, float *frac);
float NR_fmin(float x[]);
void NR_four1(float data[], unsigned long nn, int isign);
void NR_fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void NR_fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void NR_fourn(float data[], unsigned long nn[], int ndim, int isign);
void NR_fpoly(float x, float p[], int np);
void NR_fred2(int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
float NR_fredin(float x, int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
void NR_frenel(float x, float *s, float *c);
void NR_frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
void NR_ftest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *f, float *prob);
float NR_gamdev(int ia, long *idum);
float NR_gammln(float xx);
float NR_gammp(float a, float x);
float NR_gammq(float a, float x);
float NR_gasdev(long *idum);
void NR_gaucof(int n, float a[], float b[], float amu0, float x[], float w[]);
void NR_gauher(float x[], float w[], int n);
void NR_gaujac(float x[], float w[], int n, float alf, float bet);
void NR_gaulag(float x[], float w[], int n, float alf);
void NR_gauleg(float x1, float x2, float x[], float w[], int n);
void NR_gaussj(float **a, int n, float **b, int m);
void NR_gcf(float *gammcf, float a, float x, float *gln);
float NR_golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);
void NR_gser(float *gamser, float a, float x, float *gln);
void NR_hpsel(unsigned long m, unsigned long n, float arr[], float heap[]);
void NR_hpsort(unsigned long n, float ra[]);
void NR_hqr(float **a, int n, float wr[], float wi[]);
void NR_hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void NR_hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void NR_hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode);
void NR_hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode);
void NR_hunt(float xx[], unsigned long n, float x, unsigned long *jlo);
void NR_hypdrv(float s, float yy[], float dyyds[]);
nr_fcomplex NR_hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void NR_hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short NR_icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short NR_icrc1(unsigned short crc, unsigned char onech);
unsigned long NR_igray(unsigned long n, int is);
void NR_iindexx(unsigned long n, long arr[], unsigned long indx[]);
void NR_indexx(unsigned long n, float arr[], unsigned long indx[]);
void NR_interp(double **uf, double **uc, int nf);
int NR_irbit1(unsigned long *iseed);
int NR_irbit2(unsigned long *iseed);
void NR_jacobi(float **a, int n, float d[], float **v, int *nrot);
void NR_jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
long NR_julday(int mm, int id, int iyyy);
void NR_kendl1(float data1[], float data2[], unsigned long n, float *tau, float *z,
	float *prob);
void NR_kendl2(float **tab, int i, int j, float *tau, float *z, float *prob);
void NR_kermom(double w[], double y, int m);
void NR_ks2d1s(float x1[], float y1[], unsigned long n1,
	void (*quadvl)(float, float, float *, float *, float *, float *),
	float *d1, float *prob);
void NR_ks2d2s(float x1[], float y1[], unsigned long n1, float x2[], float y2[],
	unsigned long n2, float *d, float *prob);
void NR_ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob);
void NR_kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *d, float *prob);
void NR_laguer(fcomplex a[], int m, fcomplex *x, int *its);
void NR_lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	int ma, float **covar, float *chisq, void (*funcs)(float, float [], int));
void NR_linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	 int itmax, int *iter, double *err);
void NR_linmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []));
void NR_lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
	 float *f, float stpmax, int *check, float (*func)(float []));
void NR_load(float x1, float v[], float y[]);
void NR_load1(float x1, float v1[], float y[]);
void NR_load2(float x2, float v2[], float y[]);
void NR_locate(float xx[], unsigned long n, float x, unsigned long *j);
void NR_lop(double **out, double **u, int n);
void NR_lubksb(float **a, int n, int *indx, float b[]);
void NR_ludcmp(float **a, int n, int *indx, float *d);
void NR_machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	float *eps, float *epsneg, float *xmin, float *xmax);
void NR_matadd(double **a, double **b, double **c, int n);
void NR_matsub(double **a, double **b, double **c, int n);
void NR_medfit(float x[], float y[], int ndata, float *a, float *b, float *abdev);
void NR_memcof(float data[], int n, int m, float *xms, float d[]);
int NR_metrop(float de, float t);
void NR_mgfas(double **u, int n, int maxcyc);
void NR_mglin(double **u, int n, int ncycle);
float NR_midexp(float (*funk)(float), float aa, float bb, int n);
float NR_midinf(float (*funk)(float), float aa, float bb, int n);
float NR_midpnt(float (*func)(float), float a, float b, int n);
float NR_midsql(float (*funk)(float), float aa, float bb, int n);
float NR_midsqu(float (*funk)(float), float aa, float bb, int n);
void NR_miser(float (*func)(float []), float regn[], int ndim, unsigned long npts,
	float dith, float *ave, float *var);
void NR_mmid(float y[], float dydx[], int nvar, float xs, float htot,
	int nstep, float yout[], void (*derivs)(float, float[], float[]));
void NR_mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
	float *fc, float (*func)(float));
void NR_mnewt(int ntrial, float x[], int n, float tolx, float tolf);
void NR_moment(float data[], int n, float *ave, float *adev, float *sdev,
	float *var, float *skew, float *curt);
void NR_mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void NR_mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void NR_mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void NR_mpinv(unsigned char u[], unsigned char v[], int n, int m);
void NR_mplsh(unsigned char u[], int n);
void NR_mpmov(unsigned char u[], unsigned char v[], int n);
void NR_mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void NR_mpneg(unsigned char u[], int n);
void NR_mppi(int n);
void NR_mprove(float **a, float **alud, int n, int indx[], float b[],
	float x[]);
void NR_mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void NR_mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void NR_mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void NR_mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void NR_mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n);
void NR_mrqcof(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int));
void NR_mrqmin(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda);
void NR_newt(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void NR_odeint(float ystart[], int nvar, float x1, float x2,
	float eps, float h1, float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*rkqs)(float [], float [], int, float *, float, float,
	float [], float *, float *, void (*)(float, float [], float [])));
void NR_orthog(int n, float anu[], float alpha[], float beta[], float a[],
	float b[]);
void NR_pade(double cof[], int n, float *resid);
void NR_pccheb(float d[], float c[], int n);
void NR_pcshft(float a, float b, float d[], int n);
void NR_pearsn(float x[], float y[], unsigned long n, float *r, float *prob,
	float *z);
void NR_period(float x[], float y[], int n, float ofac, float hifac,
	float px[], float py[], int np, int *nout, int *jmax, float *prob);
void NR_piksr2(int n, float arr[], float brr[]);
void NR_piksrt(int n, float arr[]);
void NR_pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	float ***c, float **s);
float NR_plgndr(int l, int m, float x);
float NR_poidev(float xm, long *idum);
void NR_polcoe(float x[], float y[], int n, float cof[]);
void NR_polcof(float xa[], float ya[], int n, float cof[]);
void NR_poldiv(float u[], int n, float v[], int nv, float q[], float r[]);
void NR_polin2(float x1a[], float x2a[], float **ya, int m, int n,
	float x1, float x2, float *y, float *dy);
void NR_polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void NR_powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	float (*func)(float []));
void NR_predic(float data[], int ndata, float d[], int m, float future[], int nfut);
float NR_probks(float alam);
void NR_psdes(unsigned long *lword, unsigned long *irword);
void NR_pwt(float a[], unsigned long n, int isign);
void NR_pwtset(int n);
float NR_pythag(float a, float b);
void NR_pzextr(int iest, float xest, float yest[], float yz[], float dy[],
	int nv);
float NR_qgaus(float (*func)(float), float a, float b);
void NR_qrdcmp(float **a, int n, float *c, float *d, int *sing);
float NR_qromb(float (*func)(float), float a, float b);
float NR_qromo(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
void NR_qroot(float p[], int n, float *b, float *c, float eps);
void NR_qrsolv(float **a, int n, float c[], float d[], float b[]);
void NR_qrupdt(float **r, float **qt, int n, float u[], float v[]);
float NR_qsimp(float (*func)(float), float a, float b);
float NR_qtrap(float (*func)(float), float a, float b);
float NR_quad3d(float (*func)(float, float, float), float x1, float x2);
void NR_quadct(float x, float y, float xx[], float yy[], unsigned long nn,
	float *fa, float *fb, float *fc, float *fd);
void NR_quadmx(float **a, int n);
void NR_quadvl(float x, float y, float *fa, float *fb, float *fc, float *fd);
float NR_ran0(long *idum);
float NR_ran1(long *idum);
float NR_ran2(long *idum);
float NR_ran3(long *idum);
float NR_ran4(long *idum);
void NR_rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void NR_ranpt(float pt[], float regn[], int n);
void NR_ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
void NR_ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev);
double NR_ratval(double x, double cof[], int mm, int kk);
float NR_rc(float x, float y);
float NR_rd(float x, float y, float z);
void NR_realft(float data[], unsigned long n, int isign);
void NR_rebin(float rc, int nd, float r[], float xin[], float xi[]);
void NR_red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, float ***c, float **s);
void NR_relax(double **u, double **rhs, int n);
void NR_relax2(double **u, double **rhs, int n);
void NR_resid(double **res, double **u, double **rhs, int n);
float NR_revcst(float x[], float y[], int iorder[], int ncity, int n[]);
void NR_reverse(int iorder[], int ncity, int n[]);
float NR_rf(float x, float y, float z);
float NR_rj(float x, float y, float z, float p);
void NR_rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	void (*derivs)(float, float [], float []));
void NR_rkck(float y[], float dydx[], int n, float x, float h,
	float yout[], float yerr[], void (*derivs)(float, float [], float []));
void NR_rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []));
void NR_rkqs(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void NR_rlft3(float ***data, float **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
float NR_rofunc(float b);
void NR_rotate(float **r, float **qt, int n, int i, float a, float b);
void NR_rsolv(float **a, int n, float d[], float b[]);
void NR_rstrct(double **uc, double **uf, int nc);
float NR_rtbis(float (*func)(float), float x1, float x2, float xacc);
float NR_rtflsp(float (*func)(float), float x1, float x2, float xacc);
float NR_rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float NR_rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float NR_rtsec(float (*func)(float), float x1, float x2, float xacc);
void NR_rzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);
void NR_savgol(float c[], int np, int nl, int nr, int ld, int m);
void NR_score(float xf, float y[], float f[]);
void NR_scrsho(float (*fx)(float));
float NR_select(unsigned long k, unsigned long n, float arr[]);
float NR_selip(unsigned long k, unsigned long n, float arr[]);
void NR_shell(unsigned long n, float a[]);
void NR_shoot(int n, float v[], float f[]);
void NR_shootf(int n, float v[], float f[]);
void NR_simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
	float *bmax);
void NR_simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1);
void NR_simp3(float **a, int i1, int k1, int ip, int kp);
void NR_simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void NR_simpr(float y[], float dydx[], float dfdx[], float **dfdy,
	int n, float xs, float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void NR_sinft(float y[], int n);
void NR_slvsm2(double **u, double **rhs);
void NR_slvsml(double **u, double **rhs);
void NR_sncndn(float uu, float emmc, float *sn, float *cn, float *dn);
double NR_snrm(unsigned long n, double sx[], int itol);
void NR_sobseq(int *n, float x[]);
void NR_solvde(int itmax, float conv, float slowc, float scalv[],
	int indexv[], int ne, int nb, int m, float **y, float ***c, float **s);
void NR_sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac);
void NR_sort(unsigned long n, float arr[]);
void NR_sort2(unsigned long n, float arr[], float brr[]);
void NR_sort3(unsigned long n, float ra[], float rb[], float rc[]);
void NR_spctrm(FILE *fp, float p[], int m, int k, int ovrlap);
void NR_spear(float data1[], float data2[], unsigned long n, float *d, float *zd,
	float *probd, float *rs, float *probrs);
void NR_sphbes(int n, float x, float *sj, float *sy, float *sjp, float *syp);
void NR_splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void NR_splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y);
void NR_spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void NR_splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void NR_spread(float y, float yy[], unsigned long n, float x, int m);
void NR_sprsax(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void NR_sprsin(float **a, int n, float thresh, unsigned long nmax, float sa[],
	unsigned long ija[]);
void NR_sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float sc[], unsigned long ijc[]);
void NR_sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float thresh, unsigned long nmax, float sc[], unsigned long ijc[]);
void NR_sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[]);
void NR_sprstx(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void NR_stifbs(float y[], float dydx[], int nv, float *xx,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void NR_stiff(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void NR_stoerm(float y[], float d2y[], int nv, float xs,
	float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void NR_svbksb(float **u, float w[], float **v, int m, int n, float b[],
	float x[]);
void NR_svdcmp(float **a, int m, int n, float w[], float **v);
void NR_da_svdcmp2(float **a, int m, int n, float w[], float **v, int maxiter);
void NR_da_svdcmp3(float **a, int m, int n, float w[], float **v, int maxiter);
void NR_svdfit(float x[], float y[], float sig[], int ndata, float a[],
	int ma, float **u, float **v, float w[], float *chisq,
	void (*funcs)(float, float [], int));
void NR_svdvar(float **v, int ma, float w[], float **cvm);
void NR_toeplz(float r[], float x[], float y[], int n);
void NR_tptest(float data1[], float data2[], unsigned long n, float *t, float *prob);
void NR_tqli(float d[], float e[], int n, float **z);
float NR_trapzd(float (*func)(float), float a, float b, int n);
void NR_tred2(float **a, int n, float d[], float e[]);
void NR_tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n);
float NR_trncst(float x[], float y[], int iorder[], int ncity, int n[]);
void NR_trnspt(int iorder[], int ncity, int n[]);
void NR_ttest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void NR_tutest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void NR_twofft(float data1[], float data2[], float fft1[], float fft2[],
	unsigned long n);
void NR_vander(double x[], double w[], double q[], int n);
void NR_vegas(float regn[], int ndim, float (*fxn)(float [], float), int init,
	unsigned long ncall, int itmx, int nprn, float *tgral, float *sd,
	float *chi2a);
void NR_voltra(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(int, float), float (*ak)(int, int, float, float));
void NR_wt1(float a[], unsigned long n, int isign,
	void (*wtstep)(float [], unsigned long, int));
void NR_wtn(float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(float [], unsigned long, int));
void NR_wwghts(float wghts[], int n, float h,
	void (*kermom)(double [], double ,int));
int NR_zbrac(float (*func)(float), float *x1, float *x2);
void NR_zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb);
float NR_zbrent(float (*func)(float), float x1, float x2, float tol);
void NR_zrhqr(float a[], int m, float rtr[], float rti[]);
float NR_zriddr(float (*func)(float), float x1, float x2, float xacc);
void NR_zroots(fcomplex a[], int m, fcomplex roots[], int polish);

#endif /* _NR_H_ */
