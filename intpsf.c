#include <stdio.h>
#include <math.h>

#ifndef MACOSX
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include "tinytim.h"

#ifdef TT_THREADED
#include <pthread.h>
#endif

static float **Sinc_table;
#define K 21
#define DK 10

struct iparstruct {
	float 	mag;
        int   	size_out;
        float 	**psf_out;
        int   	size_in;
        complex **psf_in;
	int   	lower;
	int   	upper;
};

struct tparstruct {
        int start_row;
        struct iparstruct *ipars;
};

/*****************************************************************
  File : intpsf.c

  Contents :
        Pixel_mtf : convolve PSF with pixel response function
  Convolve_jitter : convolve PSF with Gaussian jitter
         Sinc_int : sinc interpolation
    Integrate_psf : integrate Nyquist sampled PSF over detector pixels
  
  These routines convolve an input PSF image with the detector
  response function for an idealized square detector pixel and
  resamples the result.  The end product is a PSF integrated over
  detector pixels at detector pixel spacings.

  The input PSF is assumed to be critically sampled (or very near
  so).  Sinc interpolation is used for sampling.

  John Krist
  January 2000 (modified for SIRTF)
******************************************************************/ 

#define SQR(x) ((x)*(x))

/*-----------------------------------------------------------------------
* Pixel_mtf :
*
*	Convolve the critically sampled PSF with the detector pixel response 
*	function.  This is equivalent to multiplying the optical transfer 
*       function (the Fourier transform of the PSF) with the detector pixel 
*       MTF in Fourier space.  Note that this assumes that the detector
*       pixel is perfectly square with uniform response over its entire area.
*
* Inputs :
*        otf :  Complex valued 2-d array containing Fourier transform of PSF.
*  psf_scale :  Scale of input PSF (units arbitrary; must be same units as
*		  "pixel_size"
*          n :  Dimension of otf array (n by n).
* pixel_size :  Pixel size in same units as psf_scale.
*    camera  :  camera id number
*
* Outputs :
*        Returns pixel-mtf-multipled otf in "otf".
*-----------------------------------------------------------------------*/
static void Pixel_mtf( complex **otf, float psf_scale, int n, float pixel_size )
{
	float	p, inc, constant, *x, t;
	int	i, j;


	x = (float *)malloc( n * sizeof(float) );
	p = 0.5 * pixel_size / psf_scale;

	/* Generate a lookup table containing the pixel MTF */

	inc = 2.0 / n;
	constant = M_PI * p;

	x[0] = 0.0;
	x[n/2] = -1.;

	for ( i = 1; i < n/2; ++i )
	{
		x[i] = x[i-1] + inc;
		x[i+n/2] = x[i-1+n/2] + inc;
	}

	for ( i = 0; i < n; ++i )
		x[i] = sinc( x[i] * constant );


	/* Convolve PSF with pixel response function     *
 	 * (multiply OTF by pixel MTF in Fourier domain) */

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < n; ++i )
		{
			t = x[i] * x[j];
			otf[j][i].r *= t;
			otf[j][i].i *= t;
		}
	}

	free( x );

} /* Pixel_mtf */

/*---------------------------------------------------------------------
*  Convolve_jitter :
*	Convolve the PSF with a jitter function.
*
*  Inputs :
*      psf_in : complex array containing OTF (FFT) of Nyquist-sampled PSF
*    scale_in : pixel scale of psf_in
*     size_in : dimensions of psf_in (size_in by size_in pixels)
*      jitter : one-sigma value of jitter in same units as scale_in (arcsec)
*			
*----------------------------------------------------------------------*/
void Convolve_jitter( complex **psf_in, float scale_in, int size_in, 
		     float jitter_major, float jitter_minor, float jitter_angle )
{
	float	cx, cy, v, cos_angle, sin_angle, radius, xj0, yj0, xj, yj;
	int	x, y;

	radius = size_in / 2;

	cos_angle = cos(-jitter_angle * M_PI / 180);
	sin_angle = sin(-jitter_angle * M_PI / 180);

	/* Multiply the OTF (FFT of PSF) by the jitter function */

	cx = M_PI * jitter_major / (2 * scale_in);
	cy = M_PI * jitter_minor / (2 * scale_in);
        Shift_complex_to_center( psf_in, size_in );

	for ( y = 0; y < size_in; ++y )
	{
		for ( x = 0; x < size_in; ++x )
		{
			xj0 = (x - radius) / radius;
			yj0 = (y - radius) / radius;

			xj = xj0*cos_angle - yj0*sin_angle;
			yj = xj0*sin_angle + yj0*cos_angle;
			xj = xj * cx;
			yj = yj * cy;

			v = exp( -2.0 * (xj*xj + yj*yj) );
			psf_in[y][x].r *= v;
			psf_in[y][x].i *= v;
		}
	}

        Shift_complex_to_center( psf_in, size_in );

} /* Convolve_jitter */

/*------------------------------------------------------------------------------------
*  Convolve_jitter_old :
*	Convolve the Nyquist-sampled PSF with a Gaussian jitter function by
*	multiplying the OTF (FFT of PSF) by the jitter MTF in Fourier space.
*
*  Inputs :
*      psf_in : complex array containing OTF (FFT) of Nyquist-sampled PSF
*    scale_in : pixel scale of psf_in
*     size_in : dimensions of psf_in (size_in by size_in pixels)
*      jitter : one-sigma value of jitter in same units as scale_in (arcsec)
------------------------------------------------------------------------------------*/
static void Convolve_jitter_old( complex **psf_in, float scale_in, int size_in, float jitter )
{
        float   constant, radius, dist, dist_sqr[4096], y_sqr, v;
        int     x, y;


        /* Create a lookup table of distances.  The table is shifted *
         * to account for the shifting of the OTF.  This essentially *
         * creates a jitter function shifted to (0,0).               */

        radius = size_in / 2;

        for ( x = 0; x < size_in/2; ++x )
        {
                dist = x / radius;
                dist_sqr[x] = dist * dist;
        }

        for ( x = size_in/2; x < size_in; ++x )
        {
                dist = (size_in - x) / radius;
                dist_sqr[x] = dist * dist;
        }

        /* Multiply the OTF (FFT of PSF) by the jitter MTF */

	constant = M_PI * jitter / (2 * scale_in);

        for ( y = 0; y < size_in; ++y )
        {
                y_sqr = dist_sqr[y];

                for ( x = 0; x < size_in; ++x )
                {
			v = sqrt(dist_sqr[x] + y_sqr) * constant;
                        v = exp( -2.0 * v * v );
                        psf_in[y][x].r *= v;
                        psf_in[y][x].i *= v;
                }
        }

} /* Convolve_jitter */

/*-----------------------------------------------------------------------------*/
static void Int_row( int row, struct iparstruct *ipars )
{
        float   x_in, y_in, sum;
        int     column, i, j, i0, j0, x, y;

        y_in = (row - ipars->size_out/2) * ipars->mag;
        j0 = (int)(y_in + ipars->size_in/2);
        if ( j0 < ipars->lower || j0 > ipars->upper )
                return;

        for ( column = 0; column < ipars->size_out; ++column )
        {
                x_in = (column - ipars->size_out/2) * ipars->mag;
                i0 = (int)(x_in + ipars->size_in/2);
                if ( i0 < ipars->lower || i0 > ipars->upper )
                        continue;

                sum = 0.0;

                for ( j = -DK; j <= DK; ++j)
                {
                        y = j0 + j;
                        for ( i = -DK; i <= DK; ++i )
                        {
                                x = i0 + i;
                                sum += ipars->psf_in[y][x].r * Sinc_table[column][i+DK] * Sinc_table[row][j+DK];
                        }
                }

		if ( sum > 0.0 )
                	ipars->psf_out[row][column] = sum;
		else
                	ipars->psf_out[row][column] = 0.0;
        }

}


#ifdef TT_THREADED
/*---------------------------------------------------------------------------
* Integrate_thread :
*---------------------------------------------------------------------------*/
static void *Integrate_thread( void *vpars )
{
        int     row;
        struct  tparstruct *tpars = vpars;

        for ( row = tpars->start_row; row < tpars->ipars->size_out; row += Num_threads )
                Int_row( row, tpars->ipars );

        return( NULL );
}
#endif

/*-----------------------------------------------------------------------------*/
static void Sinc_int( complex **psf_in, float scale_in, int size_in,
               float **psf_out, float scale_out, int size_out )
{
#ifdef TT_THREADED
        pthread_t thread[MAX_THREADS];
	struct  tparstruct tpars[MAX_THREADS];
	int	ithread;
	void	*retval;
#endif
        int     pix_out, y, i;
        float   pix_in, r;
        struct  iparstruct ipars;


        ipars.mag = scale_out / scale_in;
        ipars.size_out = size_out;
        ipars.psf_out = psf_out;
        ipars.size_in = size_in;
        ipars.psf_in = psf_in;
        ipars.lower = DK + 2;
        ipars.upper = size_in - DK - 2;

        /* create lookup table of sinc interpolator coefficients */

        Sinc_table = Alloc_image( K, size_out );

        for ( pix_out = 0; pix_out < size_out; ++pix_out )
        {
                pix_in = (pix_out - size_out/2) * ipars.mag + size_in/2;

                for ( i = -DK; i <= DK; ++i )
                {
                        r = i + (int)pix_in - pix_in;
                        if ( fabs(r) <= DK )
                                Sinc_table[pix_out][i+DK] = sinc(M_PI*r) * sinc(M_PI*r/DK);
                        else
                                Sinc_table[pix_out][i+DK] = 0.0;
                }
        }

#ifndef TT_THREADED
	for ( y = 0; y < size_out; ++y )
		Int_row( y, &ipars );
#else
        for ( ithread = 0; ithread < Num_threads; ++ithread )
        {
                tpars[ithread].start_row = ithread;
                tpars[ithread].ipars = &ipars;
                if ( pthread_create(&thread[ithread], NULL, Integrate_thread, &tpars[ithread]) )
                {
                        fprintf(stderr, "Cannot make thread %d\n", ithread);
                        exit(1);
                }
        }

        /* Join (collapse) the two threads */

        for ( ithread = 0; ithread < Num_threads; ++ithread )
        {
                if ( pthread_join(thread[ithread], &retval) )
                {
                        fprintf(stderr, "Thread join failed\n");
                        exit(1);
                }
        }
#endif

        Free_image( Sinc_table );

} /* Sinc_int */

/*----------------------------------------------------------------------------
*  Integrate_psf :
*
*     Integrate the Nyquist (critically) sampled PSF "psf_in" (with a sampling
*     of "scale_in") onto "scale_out" sized detector pixels.  This is done
*     by convolving the critically sampled PSF with the pixel response
*     function and then sampling the result at detector-pixel spacings.
*
*  Inputs :
*      psf_in : Complex valued 2D array containing critically sampled PSF;
*		   (only the real portion of this data is used);
*    scale_in : Spacing of "psf_in" values in same units as "scale_out"
*     size_in : Dimension of "psf_in" (size_in by size_in)
*     psf_out : 2D array allocated by the caller into which the final
*	  	   PSF is placed
*   scale_out : Spacing of "psf_out" values (detector pixel size) in same 
*		   units as "scale_in"
*    size_out : Dimension of "psf_out"
*jitter_major : Amount of Gaussian jitter along major jitter axis in RMS mas
*jitter_minor : Amount of Gaussian jitter along minor jitter axis in RMS mas
*jitter_angle : Direction of major jitter axis from +x in degrees
* Returns :
*    Places integrated PSF in "psf_out".
*---------------------------------------------------------------------------*/ 
void Integrate_psf( complex **psf_in, float scale_in, int size_in,
		   float **psf_out, float scale_out, int size_out,
		   float jitter_major, float jitter_minor, float jitter_angle )
{
	float	i, r;
	int	x, y;

	/* compute Fourier transform of PSF (=OTF) */

        fft2d( psf_in[0], size_in, -1 );

	/* multiply OTF with pixel MTF */

	Pixel_mtf( psf_in, scale_in, size_in, scale_out );

	/* if jitter specified, convolve PSF with jitter  *
	 * (multiply PSF by jitter MTF in Fourier space)  */

	if ( jitter_major > 0.001 )   /* jitter is given in milliarcsec */
		Convolve_jitter( psf_in, scale_in, size_in, jitter_major*1.0e-3, jitter_minor*1.0e-3, jitter_angle );
	
	/* Transform MTF back into PSF */

	fft2d( psf_in[0], size_in, 1 );

        for ( y = 0; y < size_in; ++y )
        {
                for ( x = 0; x < size_in; ++x )
                {
                        r = psf_in[y][x].r;
                        i = psf_in[y][x].i;

			/* ensure PSF is real-valued; FFT noise can often *
			 * creep into the imaginary portion               */

                        psf_in[y][x].r = sqrt(r*r + i*i);
			psf_in[y][x].i = 0.0;
                }
        }

        Shift_complex_to_center( psf_in, size_in );

	/* Sample convolved PSF using sinc interpolation */

	Sinc_int( psf_in, scale_in, size_in, psf_out, scale_out, size_out );

} /* Integrate_psf */

