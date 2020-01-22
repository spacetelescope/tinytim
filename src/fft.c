/* File     :  fft.c
 * 
 * Contents : 
 *	fft  :  Compute the FFT of a multidimensional array.
 *	fft2d : Compute FFT of a 2-d array.
 *
 * Author   :  Converted from Fortran to C by John Krist (STScI)
 *	       Original author unknown (Singleton?).  Obtained from 
 *	       netlib.
 * Date     :  January 1992
 *
 * Modifications:
 *   August 1994 - JEK
 *	Modified routines to use complex structures, which improves cache
 *	hits.  Thanks to David Robinson for the suggestion.	
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*---------------------------------------------------------------------*/
static void fft( complex *image, int ntot, int n, int nspan, int isn )
{
  /*----------------------------------------------------------------------
  Converted from Fortran to C (by brute force) by John Krist, STScI.
  This code was obtained from the GO directory on netlib.  The calls for
  it are the same as in Fortran.  The following description is from the
  original source:

  Multivariate complex fourier transform, computed in place using mixed-radix 
  fast Fourier transform algorithm.
  By R. C. Singleton, Stanford Research Institute, Sept. 1968.
  Arrays a and b originally hold the real and imaginary components of the data,
  and return the real and imaginary components of the resulting Fourier 
  coefficients.  Multivariate data is indexed according to the Fortran array 
  element successor function, without limit on the number of implied multiple 
  subscripts.  The subroutine is called once for each variate.  The calls for 
  a multivariate transform may be in any order.  

  ntot is the total number of complex data values.
  n is the dimension of the current variable.
  nspan/n is the spacing of consecutive data values while indexing the current 
	variable.
  The sign of isn determines the sign of the complex exponential, and the 
	magnitude of isn is normally one.
  A tri-variate transform with a(n1,n2,n3), b(n1,n2,n3) is computed by
      call fft(a,b,n1*n2*n3,n1,n1,1)
      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
  For a single-variate transform,
	 ntot = n = nspan = (number of complex data values), e.g.
      		call fft(a,b,n,n,n,1)
  The data can alternatively be stored in a single complex array c in standard 
  Fortran fashion, i.e. alternating real and imaginary parts. then with most 
  Fortran compilers, the complex array c can be equivalenced to a real array a,
  the magnitude of isn changed to two to give correct indexing increment, and 
  a(1) and a(2) used to pass the initial addresses for the sequences of real 
  and imaginary values, e.g.
       complex c(ntot)
       real    a(2*ntot)
       equivalence (c(1),a(1))
       call fft(a(1),a(2),ntot,n,nspan,2)
  Arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp) are used for 
  temporary storage.  if the available storage is insufficient, the program is 
  terminated by a stop.  maxf must be >= the maximum prime factor of n.
  maxp must be > the number of prime factors of n.  In addition, if the 
  square-free portion k of n has two or more prime factors, then maxp must 
  be >= k-1.
     dimension a(1),b(1)
  Array storage in nfac for a maximum of 15 prime factors of n.
  If n has more than one square-free factor, the product of the
    square-free factors must be <= 210
  -------------------------------------------------------------------------*/

#define ii i

/*  MAXPRIME is the (largest prime number + 1) that can be handled */

#define MAXPRIME  128

	int        nfac[12], np[210];
	int        i, j, jj, k, kk, k1, k2, k3, k4, m, jf, maxp, inc;
	int        nt, ks, kspan, kspnn, nn, jc, kt, maxf;
  
/* array storage for maximum prime factor of MAXPRIME-1 */

	float      at[MAXPRIME], ck[MAXPRIME], bt[MAXPRIME], sk[MAXPRIME];
	float      c72, s72, s120, rad, radf, sd, cd, ak, bk, s1, c1;
	float      aj, bj, akp, ajp, akm, bkp, bjp, bkm, c2, c3, s2, s3;
	float      ajm, bjm, bb, aa, tt;    

        c2 = c3 = s2 = s3 = 0.0;

/* This program indexes everything according to Fortran notation (first array
 * element is index=1.  Decrement the data pointer so that this will work.  */
	--image;

/* the following two constants should agree with the array dimensions. */

	k3 = 0;
	maxf = MAXPRIME - 1;
	maxp = 127;
	if (n < 2) return;
	inc = isn;
	c72 = 0.30901699437494742;
	s72 = 0.95105651629515357;
	s120 = 0.86602540378443865;
	rad = 6.2831853071796;
	if (isn < 0) 
	{
		s72 = -s72;
		s120 = -s120;
		rad = -rad;
		inc = -inc;
	}
	nt = inc * ntot;
	ks = inc * nspan;
	kspan = ks;
	nn = nt - inc;
	jc = ks / n;
	radf = rad * (float)jc * 0.5;
	i = 0;
	jf = 0;

/*  determine the factors of n */

	m = 0;
	k = n;

	while ( k-(k/16)*16 == 0 )
	{
		m = m + 1;
		nfac[m] = 4;
		k = k / 16;
	}

	j = 3;
	jj = 9;

	do {
		while ( k % jj == 0 )
		{
			m = m + 1;
			nfac[m] = j;
			k = k / jj;
		}
		j = j + 2;
		jj = j * j;
	} while (jj <= k); 

	if (k <= 4) 
	{
		kt = m;
		nfac[m+1] = k;
		if (k != 1) m = m + 1;
	}
	else
	{
		if (k-(k/4)*4 == 0)
		{
			m = m + 1;
			nfac[m] = 2;
			k = k / 4;
		}

		kt = m;
		j = 2;

                do {
			if (k % j == 0)
			{
				m = m + 1;
				nfac[m] = j;
				k = k / j;
			}
			j = ((j+1) / 2) * 2 + 1;
		} while (j <= k);
	}

	if (kt != 0)
	{
		j = kt;
		do {
			m = m + 1;
			nfac[m] = nfac[j];
			j = j - 1;
		} while (j != 0);
	}

/*  compute fourier transform */

  L100:
	sd = radf / (float)kspan;
	tt = sin(sd);
	cd = 2.0 * tt * tt;
	sd = sin(sd+sd);
	kk = 1;
	i = i + 1;
	if (nfac[i] != 2) goto L400;

/*  transform for factor of 2 (including rotation factor) */

	kspan = kspan / 2;
	k1 = kspan + 2;

	do {
		do {
			k2 = kk + kspan;
			ak = image[k2].r;
			bk = image[k2].i;
			image[k2].r = image[kk].r - ak;
			image[k2].i = image[kk].i - bk;
			image[kk].r = image[kk].r + ak;
			image[kk].i = image[kk].i + bk;
			kk = k2 + kspan;
		} while (kk <= nn);
		kk = kk - nn;
	} while (kk <= jc);

	if (kk > kspan) goto L800;

	do {
	     c1 = 1.0 - cd;
	     s1 = sd;
	     do {
		   do {
			 do {
				k2 = kk + kspan;
				ak = image[kk].r - image[k2].r;
				bk = image[kk].i - image[k2].i;
				image[kk].r = image[kk].r + image[k2].r;
				image[kk].i = image[kk].i + image[k2].i;
				image[k2].r = c1 * ak - s1 * bk;
				image[k2].i = s1 * ak + c1 * bk;
				kk = k2 + kspan;
			 } while (kk < nt);
			 k2 = kk - nt;
			 c1 = -c1;
			 kk = k1 - k2;
		   } while (kk > k2);
		   ak = c1 - (cd * c1 + sd * s1);
		   s1 = (sd * c1 - cd * s1) + s1;
		   c1 = 2.0 - (ak*ak + s1*s1);
		   s1 = c1 * s1;
		   c1 = c1 * ak;
		   kk = kk + jc;
	     } while (kk < k2);
	     k1 = k1 + 2 * inc;
	     kk = (k1 - kspan) / 2 + jc;
	} while (kk <= 2*jc);

	goto L100;

/*  transform for factor of 3 (optional code) */

  L320:
	do {
		do {
			k1 = kk + kspan;
			k2 = k1 + kspan;
			ak = image[kk].r;
			bk = image[kk].i;
			aj = image[k1].r + image[k2].r;
			bj = image[k1].i + image[k2].i;
			image[kk].r = ak + aj;
			image[kk].i = bk + bj;
			ak = -0.5 * aj + ak;
			bk = -0.5 * bj + bk;
			aj = (image[k1].r - image[k2].r) * s120;
			bj = (image[k1].i - image[k2].i) * s120;
			image[k1].r = ak - bj;         
			image[k1].i = bk + aj;
			image[k2].r = ak + bj;
			image[k2].i = bk - aj;
			kk = k2 + kspan;
		} while (kk < nn);
		kk = kk - nn;
	} while (kk <= kspan);
	
	goto L700;

/*  transform for factor of 4 */

  L400: 
      if (nfac[i] != 4) goto L600;
      kspnn=kspan;
      kspan=kspan/4;

  L410: 
      c1=1.0;
      s1=0;

  L420: 
      k1=kk+kspan;
      k2=k1+kspan;
      k3=k2+kspan;
      akp=image[kk].r+image[k2].r;
      akm=image[kk].r-image[k2].r;
      ajp=image[k1].r+image[k3].r;
      ajm=image[k1].r-image[k3].r;
      image[kk].r=akp+ajp;
      ajp=akp-ajp;
      bkp=image[kk].i+image[k2].i;
      bkm=image[kk].i-image[k2].i;
      bjp=image[k1].i+image[k3].i;
      bjm=image[k1].i-image[k3].i;
      image[kk].i=bkp+bjp;
      bjp=bkp-bjp;
      if (isn < 0) goto L450;
      akp=akm-bjm;
      akm=akm+bjm;
      bkp=bkm+ajm;
      bkm=bkm-ajm;
      if (s1 == 0) goto L460;

  L430: 
      image[k1].r=akp*c1-bkp*s1;
      image[k1].i=akp*s1+bkp*c1;
      image[k2].r=ajp*c2-bjp*s2;
      image[k2].i=ajp*s2+bjp*c2;
      image[k3].r=akm*c3-bkm*s3;
      image[k3].i=akm*s3+bkm*c3;
      kk=k3+kspan;
      if (kk <= nt) goto L420;

  L440: 
      c2=c1-(cd*c1+sd*s1);
      s1=(sd*c1-cd*s1)+s1;
      c1=2.0-(c2*c2 + s1*s1);
      s1=c1*s1;
      c1=c1*c2;
      c2=c1*c1-s1*s1;
      s2=2.0*c1*s1;
      c3=c2*c1-s2*s1;
      s3=c2*s1+s2*c1;
      kk=kk-nt+jc;
      if (kk <= kspan) goto L420;
      kk=kk-kspan+inc;
      if (kk <= jc) goto L410;
      if (kspan == jc) goto L800;
      goto L100;

  L450: 
      akp=akm+bjm;
      akm=akm-bjm;
      bkp=bkm-ajm;
      bkm=bkm+ajm;
      if (s1 != 0) goto L430;

  L460: 
      image[k1].r=akp;
      image[k1].i=bkp;
      image[k2].r=ajp;
      image[k2].i=bjp;
      image[k3].r=akm;
      image[k3].i=bkm;
      kk=k3+kspan;
      if (kk <= nt) goto L420;
      goto L440;

/*  transform for factor of 5 (optional code) */

  L510: 
      c2=c72*c72-s72*s72;
      s2=2.0*c72*s72;

  L520: 
      k1=kk+kspan;
      k2=k1+kspan;
      k3=k2+kspan;
      k4=k3+kspan;
      akp=image[k1].r+image[k4].r;
      akm=image[k1].r-image[k4].r;
      bkp=image[k1].i+image[k4].i;
      bkm=image[k1].i-image[k4].i;
      ajp=image[k2].r+image[k3].r;
      ajm=image[k2].r-image[k3].r;
      bjp=image[k2].i+image[k3].i;
      bjm=image[k2].i-image[k3].i;
      aa=image[kk].r;
      bb=image[kk].i;
      image[kk].r=aa+akp+ajp;
      image[kk].i=bb+bkp+bjp;
      ak=akp*c72+ajp*c2+aa;
      bk=bkp*c72+bjp*c2+bb;
      aj=akm*s72+ajm*s2;
      bj=bkm*s72+bjm*s2;
      image[k1].r=ak-bj;
      image[k1].i=bk+aj;
      image[k4].r=ak+bj;
      image[k4].i=bk-aj;
      ak=akp*c2+ajp*c72+aa;
      bk=bkp*c2+bjp*c72+bb;
      aj=akm*s2-ajm*s72;
      bj=bkm*s2-bjm*s72;
      image[k2].r=ak-bj;
      image[k2].i=bk+aj;
      image[k3].r=ak+bj;
      image[k3].i=bk-aj;
      kk=k4+kspan;
      if (kk < nn) goto L520;
      kk=kk-nn;
      if (kk <= kspan) goto L520;
      goto L700;

/*  transform for odd factors */

  L600: 
      k=nfac[i];
      kspnn=kspan;
      kspan=kspan/k;
      if (k == 3) goto L320;
      if (k == 5) goto L510;
      if (k == jf) goto L640;
      jf=k;
      s1=rad / (float)k;
      c1=cos(s1);
      s1=sin(s1);
      if (jf > maxf) goto L998;
      ck[jf]=1.0;
      sk[jf]=0.0;
      j=1;

	do {
		ck[j]=ck[k]*c1+sk[k]*s1;
		sk[j]=ck[k]*s1-sk[k]*c1;
		k=k-1;
		ck[k]=ck[j];
		sk[k] = -sk[j];
		j=j+1;
	} while (j < k); 

  L640: 
      k1=kk;
      k2=kk+kspnn;
      aa=image[kk].r;
      bb=image[kk].i;
      ak=aa;
      bk=bb;
      j=1;
      k1=k1+kspan;

	do {
		k2=k2-kspan;
		j=j+1;
		at[j]=image[k1].r+image[k2].r;
		ak=at[j]+ak;
		bt[j]=image[k1].i+image[k2].i;
		bk=bt[j]+bk;
		j=j+1;
		at[j]=image[k1].r-image[k2].r;
		bt[j]=image[k1].i-image[k2].i;
		k1=k1+kspan;
	} while (k1 < k2); 

      image[kk].r=ak;
      image[kk].i=bk;
      k1=kk;
      k2=kk+kspnn;
      j=1;

  L660: 
      k1=k1+kspan;
      k2=k2-kspan;
      jj=j;
      ak=aa;
      bk=bb;
      aj=0.0;
      bj=0.0;
      k=1;

  L670: 
      k=k+1;
      ak=at[k]*ck[jj]+ak;
      bk=bt[k]*ck[jj]+bk;
      k=k+1;
      aj=at[k]*sk[jj]+aj;
      bj=bt[k]*sk[jj]+bj;
      jj=jj+j;
      if (jj > jf) jj=jj-jf;
      if (k < jf) goto L670;
      k=jf-j;
      image[k1].r=ak-bj;
      image[k1].i=bk+aj;
      image[k2].r=ak+bj;
      image[k2].i=bk-aj;
      j=j+1;
      if (j < k) goto L660;
      kk=kk+kspnn;
      if (kk <= nn) goto L640;
      kk=kk-nn;
      if (kk <= kspan) goto L640;

/*  multiply by rotation factor (except for factors of 2 and 4) */

  L700: 
	if (i == m) goto L800;
	kk=jc+1;

	do {
	      c2=1.0-cd;
	      s1=sd;
	      do {
		    c1=c2;
		    s2=s1;
		    kk=kk+kspan;

		    do {
			do {
			      ak=image[kk].r;
			      image[kk].r=c2*ak-s2*image[kk].i;
			      image[kk].i=s2*ak+c2*image[kk].i;
			      kk=kk+kspnn;
			} while (kk <= nt);

			ak=s1*s2;
			s2=s1*c2+c1*s2;
			c2=c1*c2-ak;
			kk=kk-nt+kspan;
		    } while (kk <= kspnn);

		    c2=c1-(cd*c1+sd*s1);
		    s1=s1+(sd*c1-cd*s1);
		    c1=2.0-(c2*c2+s1*s1);
		    s1=c1*s1;
		    c2=c1*c2;
		    kk=kk-kspnn+jc;
	      } while (kk <= kspan);

	      kk=kk-kspan+jc+inc;
	} while (kk <= 2*jc);

	goto L100;

/*  permute the results to normal order---done in two stages  *
 *  permutation for square factors of n                       */

  L800: 
      np[1]=ks;
      if (kt == 0) goto L890;
      k=kt+kt+1;
      if (m < k) k=k-1;
      j=1;
      np[k+1]=jc;

	do {
		np[j+1]=np[j]/nfac[j];
		np[k]=np[k+1]*nfac[j];
		j=j+1;
		k=k-1;
	} while (j < k);

      k3=np[k+1];
      kspan=np[2];
      kk=jc+1;
      k2=kspan+1;
      j=1;
      if (n != ntot) goto L850;

/*  permutation for single-variate transform (optional code) */

  L820: 
      ak = image[kk].r;
      bk = image[kk].i;
      image[kk].r = image[k2].r;
      image[kk].i = image[k2].i;
      image[k2].r = ak;
      image[k2].i = bk;
      kk=kk+inc;
      k2=kspan+k2;
      if (k2 < ks) goto L820;

  L830: 
      k2=k2-np[j];
      j=j+1;
      k2=np[j+1]+k2;
      if (k2 > np[j]) goto L830;
      j=1;

  L840: 
      if (kk < k2) goto L820;
      kk=kk+inc;
      k2=kspan+k2;
      if (k2 < ks) goto L840;
      if (kk < ks) goto L830;
      jc=k3;
      goto L890;

/*  permutation for multivariate transform */

  L850: 
      k=kk+jc;

  L860: 
      ak = image[kk].r;
      bk = image[kk].i;
      image[kk].r = image[k2].r;
      image[kk].i = image[k2].i;
      image[k2].r = ak;
      image[k2].i = bk;
      kk=kk+inc;
      k2=k2+inc;
      if (kk < k) goto L860;
      kk=kk+ks-jc;
      k2=k2+ks-jc;
      if (kk < nt) goto L850;
      k2=k2-nt+kspan;
      kk=kk-nt+jc;
      if (k2 < ks) goto L850;

  L870: 
      k2=k2-np[j];
      j=j+1;
      k2=np[j+1]+k2;
      if (k2 > np[j]) goto L870;
      j=1;

  L880: 
      if (kk < k2) goto L850;
      kk=kk+jc;
      k2=kspan+k2;
      if (k2 < ks) goto L880;
      if (kk < ks) goto L870;
      jc=k3;

  L890: 
      if (2*kt+1 >= m) return;
      kspnn=np[kt+1];

/*  permutation for square-free factors of n */

      j=m-kt;
      nfac[j+1]=1;

  L900: 
      nfac[j]=nfac[j]*nfac[j+1];
      j=j-1;
      if (j != kt) goto L900;
      kt=kt+1;
      nn=nfac[kt]-1;
      if (nn > maxp) goto L998;
      jj=0;
      j=0;
      goto L906;

  L902: 
      jj=jj-k2;
      k2=kk;
      k=k+1;
      kk=nfac[k];

  L904: 
      jj=kk+jj;
      if (jj >= k2) goto L902;
      np[j]=jj;

  L906: 
      k2=nfac[kt];
      k=kt+1;
      kk=nfac[k];
      j=j+1;
      if (j <= nn) goto L904;

/*  determine the permutation cycles of length greater than 1 */

      j=0;
      goto L914;

  L910: 
      k=kk;
      kk=np[k];
      np[k] = -kk;
      if (kk != j) goto L910;
      k3=kk;

  L914: 
      j=j+1;
      kk=np[j];
      if (kk < 0) goto L914;
      if (kk != j) goto L910;
      np[j] = -j;
      if (j != nn) goto L914;
      maxf=inc*maxf;

/*  reorder a and b, following the permutation cycles */

      goto L950;

  L924: 
      j=j-1;
      if (np[j] < 0) goto L924;
      jj=jc;

  L926: 
      kspan=jj;
      if (jj > maxf) kspan=maxf;
      jj=jj-kspan;
      k=np[j];
      kk=jc*k+ii+jj;
      k1=kk+kspan;
      k2=0;

  L928: 
      k2=k2+1;
      at[k2]=image[k1].r;
      bt[k2]=image[k1].i;
      k1=k1-inc;
      if (k1 != kk) goto L928;

  L932: 
      k1=kk+kspan;
      k2=k1-jc*(k+np[k]);
      k = -np[k];

  L936: 
      image[k1].r = image[k2].r;
      image[k1].i = image[k2].i;
      k1=k1-inc;
      k2=k2-inc;
      if (k1 != kk) goto L936;
      kk=k2;
      if (k != j) goto L932;
      k1=kk+kspan;
      k2=0;

  L940: 
      k2=k2+1;
      image[k1].r = at[k2];
      image[k1].i = bt[k2];
      k1=k1-inc;
      if (k1 != kk) goto L940;
      if (jj != 0) goto L926;
      if (j != 1) goto L924;
  L950: 
      j=k3+1;
      nt=nt-kspnn;
      ii=nt-inc+1;
      if (nt >= 0) goto L924;
      return;

/*  error finish, insufficient array storage */

  L998: 
      isn=0;
      printf("FFT Error\n");
      return;

} /* fft */

/*--------------------------------------------------------------------
*  fft2d :
*	Compute the fft of a 2-d array.
*
*  Inputs :
*	image     : Complex valued image.
*	dim       : Dimension of "image" (Dim x Dim).
*	direction : 1 for inverse, -1 for forward.
*
*  Outputs :
*	The data in "image" is replaced by the fft.
*----------------------------------------------------------------------*/
void fft2d( complex *image, int dim, int direction )
{

	fft( image, dim*dim, dim, dim, direction );
	fft( image, dim*dim, dim, dim*dim, direction );

} /* fft2d */

