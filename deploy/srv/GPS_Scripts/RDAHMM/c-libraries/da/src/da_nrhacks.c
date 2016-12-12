/*******************************************************************************
MODULE NAME
da_nrhacks

ONE-LINE SYNOPSIS
Hacks of various functions in the Numerical Recipes package.

SCOPE OF THIS MODULE
All functions that are direction modifications of Numerical Recipes functions
should go here.

SEE ALSO
There should be no overlap with other modules; for the original functions 
reference the nr library.

REFERENCE(S)
-

PROGRAM EXAMPLE(S)
-

NOTES
-

RG
*******************************************************************************/
#ifndef lint
static char rcsid[] = "$Id: da_nrhacks.c,v 1.3 1997/08/24 00:53:18 granat Exp granat $";
#endif
/* This string variable allows the RCS identification info to be printed. */

/* 
 * $Log: da_nrhacks.c,v $
 * Revision 1.3  1997/08/24 00:53:18  granat
 * changed function names to DA_ convention
 *
 * Revision 1.2  1997/06/02 15:37:46  granat
 * changed to use new NR naming convention
 *
 * Revision 1.1  1997/05/13 20:57:55  granat
 * Initial revision
 *
 * */

/* C library */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <signal.h>

/* UT library */
#include "ut_error.h"
#include "ut_output.h"
#include "ut_string.h"
#include "ut_types.h"

/* NR library */
#include "nr.h"

/* this module's header */
#include "da_nrhacks.h"

/*******************************************************************************
DA_NRERROR
This version of the standard Numerical Recipes error handler raises a SIGSTOP
instead of exiting.
*******************************************************************************/
void DA_nrerror(char error_text[])
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now raising SIGSTOP flag...\n");
        raise(SIGSTOP);
}

/*******************************************************************************
DA_SVDCMP
Hacked version of the Numerical Recipes svdcmp function.
The maximum number of iterations is now an input parameter.
*******************************************************************************/
int DA_svdcmp(float **a, int m, int n, float w[], float **v, int maxiter)
{
  int flag,i,its,j,jj,k,l,nm;
  float anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  char tempstr[5];
 
  rv1=NR_vector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++)
        scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -NR_SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=i;k<=m;k++)
            s += a[k][i]*a[k][j];
          f=s/h;
          for (k=i;k<=m;k++)
            a[k][j] += f*a[k][i];
        }
        for (k=i;k<=m;k++)
          a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++)
        scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -NR_SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++)
          rv1[k]=a[i][k]/h;
        for (j=l;j<=m;j++) {
          for (s=0.0,k=l;k<=n;k++)
            s += a[j][k]*a[i][k];
          for (k=l;k<=n;k++)
            a[j][k] += s*rv1[k];
        }
        for (k=l;k<=n;k++)
          a[i][k] *= scale;
      }
    }
    anorm=NR_FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++)
            s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++)
            v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++)
        v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=NR_IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++)
      a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
        for (s=0.0,k=l;k<=m;k++)
          s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
        for (k=i;k<=m;k++)
          a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++)
        a[j][i] *= g;
    } else for (j=i;j<=m;j++)
      a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=maxiter;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
        nm=l-1;
        if ((float)(fabs(rv1[l])+anorm) == anorm) {
          flag=0;
          break;
        }
        if ((float)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if ((float)(fabs(f)+anorm) == anorm) break;
          g=w[i];
          h=NR_pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=1;j<=m;j++) {
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=1;j<=n;j++)
            v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == maxiter) {
  sprintf(tempstr, "no convergence in %d svdcmp iterations", maxiter);
  NR_error(tempstr);
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=NR_pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+NR_SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=NR_pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g = g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=NR_pythag(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=1;jj<=m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  NR_free_vector(rv1,1,n);

  return( UT_OK );
}


/*******************************************************************************
DA_DSVDCMP
Hacked version of the Numerical Recipes dsvdcmp function.
The maximum number of iterations is now an input parameter.
*******************************************************************************/
int DA_dsvdcmp(double **a, int m, int n, double w[], double **v, int maxiter)
{
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  char tempstr[5];
 
  rv1=NR_dvector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++)
        scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -NR_SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=i;k<=m;k++)
            s += a[k][i]*a[k][j];
          f=s/h;
          for (k=i;k<=m;k++)
            a[k][j] += f*a[k][i];
        }
        for (k=i;k<=m;k++)
          a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++)
        scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -NR_SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++)
          rv1[k]=a[i][k]/h;
        for (j=l;j<=m;j++) {
          for (s=0.0,k=l;k<=n;k++)
            s += a[j][k]*a[i][k];
          for (k=l;k<=n;k++)
            a[j][k] += s*rv1[k];
        }
        for (k=l;k<=n;k++)
          a[i][k] *= scale;
      }
    }
    anorm=NR_DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++)
            s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++)
            v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++)
        v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=NR_IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++)
      a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
        for (s=0.0,k=l;k<=m;k++)
          s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
        for (k=i;k<=m;k++)
          a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++)
        a[j][i] *= g;
    } else for (j=i;j<=m;j++)
      a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=maxiter;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
        nm=l-1;
        if ((double)(fabs(rv1[l])+anorm) == anorm) {
          flag=0;
          break;
        }
        if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if ((double)(fabs(f)+anorm) == anorm) break;
          g=w[i];
          h=NR_dpythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=1;j<=m;j++) {
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=1;j<=n;j++)
            v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == maxiter) {
  sprintf(tempstr, "no convergence in %d svdcmp iterations", maxiter);
  NR_error(tempstr);
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=NR_dpythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+NR_SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=NR_pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g = g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=NR_dpythag(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=1;jj<=m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  NR_free_dvector(rv1,1,n);

  return( UT_OK );
}


/*******************************************************************************
DA_AMOTRY
Hacked version of amotry designed to handle an input function (the function
being optimized) that has an extra parameter.
RG,CD
*******************************************************************************/
float DA_amotry( float **p, float *y, float *psum, int ndim, AMOEBA_PF func,
                 void *arg, int nargs, int ihi, int *nfunc, float fac )
{
    int         j;
    float       fac1, fac2;
    float       ytry;
    float       *ptry;
 
    ptry = NR_vector( 1, ndim );
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for( j = 1; j <= ndim; ++j )
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    ytry = func( arg, nargs, ptry, ndim );     /* Call the function. */
    ++*nfunc;
    if( ytry < y[ihi] ) {
        y[ihi] = ytry;
        for( j = 1; j <= ndim; ++j ) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
 
    NR_free_vector( ptry, 1, ndim );
    return( ytry );
}

/*******************************************************************************
DA_AMOEBA
Hacked version of amoeba designed to handle an input function (the function
being optimized) that has an extra parameter.  Also, the alpha, beta, and
gamma parameters of the simplex algorithm are defined as macros.  The
function also returns an error if the value of the cost function goes
below a certain threshold.
*******************************************************************************/
int DA_amoeba( float **p, float *y, int ndim, float ftol, AMOEBA_PF func,
               void *arg, int nargs, int *nfunc )
{
    int         inhi, ilo, ihi;
    int         i, j;
    int         mpts;
    float       ytry, ysave, sum, rtol;
    float       *psum;
 
    psum = NR_vector( 1, ndim );
    mpts = ndim + 1;
    *nfunc = 0;
    GET_PSUM;
    for(;;) {
        ilo = 1;
        if (y[1] > y[2])  {
            ihi  = 1;
            inhi = 2;
        } else {
            ihi  = 2;
            inhi = 1;
        };
 
        for( i = 1; i <= mpts; ++i ) {
            if( y[i] < y[ilo] )
                ilo = i;
            if( y[i] > y[ihi] ) {
                inhi = ihi;
                ihi = i;
            } else if ( y[i] > y[inhi] && i != ihi ) {
                inhi = i;
            }
        }
 
        rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs( y[ihi] ) + fabs(y[ilo]));
 
        if( rtol < ftol )
            break;
        if( *nfunc >= NMAX ) {
            NR_free_vector( psum, 1, ndim );
            return( 0 );
        }
        ytry = DA_amotry( p, y, psum, ndim, func, arg, nargs, ihi, nfunc, -ALPHA );
        if( ytry < y[ilo] )
            ytry = DA_amotry( p, y, psum, ndim, func, arg, nargs, ihi, nfunc, GAMMA );
        else if( ytry >= y[inhi] ) {
            ysave = y[ihi];
            ytry = DA_amotry( p, y, psum, ndim, func, arg, nargs, ihi, nfunc, BETA );
            if (ytry < COST_THRESH) {
              NR_free_vector( psum, 1, ndim );
              return( -1 );
            }
            if( ytry >= ysave ) {
                for( i = 1; i <= mpts; ++i ) {
                    if( i != ilo ) {
                        for( j = 1; j <= ndim; ++j ) {
                            psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                            p[i][j] = psum[j];
                        }
                        y[i] = func( arg, nargs, psum, ndim );
                        if (y[i] < COST_THRESH) {
                                NR_free_vector( psum, 1, ndim );
                                return( -1 );
                        }
                    }
                }
                *nfunc += ndim;
                GET_PSUM;
            }
        }
    }
    NR_free_vector( psum, 1, ndim );
    return( 1 );
}

