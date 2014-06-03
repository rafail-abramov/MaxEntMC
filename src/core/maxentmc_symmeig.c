/** This file is part of MaxEntMC, a maximum entropy algorithm with moment constraints. **/
/** Copyright (C) 2014 Rafail V. Abramov.                                               **/
/**                                                                                     **/
/** The routines in this file are adapted from "Numerical Recipes in C" by Press,       **/
/** Teukolsky, Vetterling and Flannery.                                                 **/
/**                                                                                     **/
/** This program is free software: you can redistribute it and/or modify it under the   **/
/** terms of the GNU General Public License as published by the Free Software           **/
/** Foundation, either version 3 of the License, or (at your option) any later version. **/
/**                                                                                     **/
/** This program is distributed in the hope that it will be useful, but WITHOUT ANY     **/
/** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     **/
/** PARTICULAR PURPOSE.  See the GNU General Public License for more details.           **/
/**                                                                                     **/
/** You should have received a copy of the GNU General Public License along with this   **/
/** program.  If not, see <http://www.gnu.org/licenses/>.                               **/

#include <math.h>
#include "maxentmc_defs.h"
#include "maxentmc_symmeig.h"

static inline void tred2(size_t const n, maxentmc_float_t * const a, size_t const tda,
                         maxentmc_float_t * const d, maxentmc_float_t * const e)
{

    size_t k,j,i;
    maxentmc_float_t scale,hh,h,g,f;

    for (i=n-1;i>0;--i) {

        h=scale=0.0;

        if (i > 1) {

            for (k=0;k<i;++k)
                scale += fabs(a[i*tda+k]);

            if (scale == 0.0)
                e[i]=a[i*tda+i-1];

            else {

                for (k=0;k<i;++k) {

                    a[i*tda+k] /= scale;
                    h += a[i*tda+k]*a[i*tda+k];
                }

                f=a[i*tda+i-1];
                g=((f >= 0.0) ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                a[i*tda+i-1]=f-g;
                f=0.0;

                for (j=0;j<i;++j) {

                    a[j*tda+i]=a[i*tda+j]/h;
                    g=0.0;

                    for (k=0;k<=j;++k)
                        g += a[j*tda+k]*a[i*tda+k];

                    for (k=j+1;k<i;++k)
                        g += a[k*tda+j]*a[i*tda+k];

                    e[j]=g/h;

                    f += e[j]*a[i*tda+j];

                }

                hh=f/(h+h);

                for (j=0;j<i;++j) {

                    f=a[i*tda+j];
                    e[j]=g=e[j]-hh*f;

                    for (k=0;k<=j;++k)
                        a[j*tda+k] -= (f*e[k]+g*a[i*tda+k]);

                }

            }

        }
        else
            e[i]=a[i*tda+i-1];

        d[i]=h;

    }

    d[0]=0.0;
    e[0]=0.0;


    for (i=0;i<n;++i){

        if (d[i]){

            for (j=0;j<i;++j){

                g=0.0;

                for (k=0;k<i;++k)
                    g += a[i*tda+k]*a[k*tda+j];

                for (k=0;k<i;++k)
                    a[k*tda+j] -= g*a[k*tda+i];

            }

        }

        d[i]=a[i*tda+i];
        a[i*tda+i]=1.0;

        for (j=0;j<i;++j)
            a[j*tda+i]=a[i*tda+j]=0.0;

    }

}

static inline maxentmc_float_t pythag(maxentmc_float_t const a, maxentmc_float_t const b)
{

    maxentmc_float_t const absa = fabs(a);
    maxentmc_float_t const absb = fabs(b);
    maxentmc_float_t temp;

    if (absa > absb){

        temp = absb/absa;
        return (absa*sqrt(1.0+temp*temp));
    }
    else{
        temp = absa/absb;
        return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+temp*temp));
    }

}

static inline int tqli(size_t const n, maxentmc_float_t * const a, size_t const tda,
                       maxentmc_float_t * const d,
                       maxentmc_float_t * const e)
{

    size_t m,l,iter,i,k;
    maxentmc_float_t s,r,p,g,f,dd,c,b;

    for (i=1;i<n;++i)
        e[i-1]=e[i];

    e[n-1]=0.0;

    for (l=0;l<n;++l) {

        iter=0;

        do {

            for (m=l;m<(n-1);++m) {

                dd=fabs(d[m])+fabs(d[m+1]);
                if ((fabs(e[m])+dd) == dd)
                    break;

            }

            if (m != l) {

                if ((++iter) == 30)
                    return -1;

                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+((g) >= 0.0 ? fabs(r) : -fabs(r)));
                s=c=1.0;
                p=0.0;

                for (i=m;i>l;--i) {

                    f=s*e[i-1];
                    b=c*e[i-1];

                    e[i]=(r=pythag(f,g));

                    if (r == 0.0) {
                        d[i] -= p;
                        e[m]=0.0;
                        break;
                    }

                    s=f/r;
                    c=g/r;
                    g=d[i]-p;
                    r=(d[i-1]-g)*s+2.0*c*b;
                    d[i]=g+(p=s*r);
                    g=c*r-b;

                    for (k=0;k<n;++k) {

                        f=a[k*tda+i];
                        a[k*tda+i]=s*a[k*tda+i-1]+c*f;
                        a[k*tda+i-1]=c*a[k*tda+i-1]-s*f;

                    }

                }

                if (r == 0.0 && i > l)
                    continue;

                d[l] -= p;
                e[l]=g;
                e[m]=0.0;

            }

        } while (m != l);

    }

    return 0;

}


static inline void eigsrt(size_t const n, maxentmc_float_t * const v,
                            size_t const tda, maxentmc_float_t * const d)
{

    size_t k,j,i;
    maxentmc_float_t p;

    if(n>1)

        for (i=0;i<(n-1);++i) {

            p=d[k=i];

            for (j=i+1;j<n;++j)
                if (d[j] > p) p=d[k=j];
                /*if (fabs(d[j]) > fabs(p)) p=d[k=j];*/

            if (k != i) {

                d[k]=d[i];
                d[i]=p;

                for (j=0;j<n;++j) {

                    p=v[j*tda+i];
                    v[j*tda+i]=v[j*tda+k];
                    v[j*tda+k]=p;

                }

            }

        }

}


int maxentmc_symmeig(size_t const n, maxentmc_float_t * const A, size_t const tda, maxentmc_float_t * const eigval)
{

    if(n==0 || tda==0){
        MAXENTMC_MESSAGE(stderr,"error: zero matrix size");
        return -1;
    }

    MAXENTMC_CHECK_NULL(A);
    MAXENTMC_CHECK_NULL(eigval);
    maxentmc_float_t temp[n];

    tred2(n, A, tda, eigval, temp);

    int status = tqli(n, A, tda, eigval, temp);

    if(status){
        MAXENTMC_MESSAGE(stderr,"error: too many iterations in tqli");
    }
    else
        eigsrt(n, A, tda, eigval);

    return status;

}
