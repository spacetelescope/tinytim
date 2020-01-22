/* File :  convotf.c
 *
 * Contents :
 *	Convolve_otf :  Convolve the critically sampled PSF with
 *		the integrated pixel OTF.
 *
 * Author   :  John Krist
 * Date     :  May 1992
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*-----------------------------------------------------------------------
* Convolve_otf :
*
*	Convolve the critically sampled PSF with the integrated
*	pixel OTF.
*
* Inputs :
*    crit_psf  :  Complex valued 2-d array containing critically sampled PSF.
*        dim  :  dimension of Real and Imaginary arrays (dim by dim).
*  data_scale  :  Pixel scale of critically sampled PSF (arcsec).
* pixel_scale  :  Pixel scale of integrated PSF (arcsec).
*
*-----------------------------------------------------------------------*/
void Convolve_otf( complex **crit_psf, int dim, float data_scale, float pixel_scale )
{
	float	pixel_size, inc, constant, x[2048], t;
	int	i, j;


	pixel_size = 0.5 * pixel_scale / data_scale;

	/* Generate a lookup table containing the integrated *
	 * pixel OTF components.			     */

	inc = 2.0 / (float)dim;
	constant = M_PI * pixel_size;

	x[0] = 0.0;
	x[dim/2] = -1.;

	for ( i = 1; i < dim/2; ++i )
	{
		x[i] = x[i-1] + inc;
		x[i+dim/2] = x[i-1+dim/2] + inc;
	}

	for ( i = 0; i < dim; ++i )
		x[i] = sinc( x[i] * constant );


	/* Convolve (multiply in Fourier space) */

	for ( j = 0; j < dim; ++j )
	{
		for ( i = 0; i < dim; ++i )
		{
			t = x[i] * x[j];
			crit_psf[j][i].r *= t;
			crit_psf[j][i].i *= t;
		}
	}

} /* Convolve_otf */

