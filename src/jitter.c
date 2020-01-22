/* File     : jitter.c
 *
 * Contents :
 *	Compute_jitter : Convolve PSF with a jitter function
 *
 * Author   : John Krist (STScI)
 * Date     : January 1992
 *
 * Modifications:
 *	May 2004 : Added 2-D elliptical Gaussian jitter
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*---------------------------------------------------------------------
*  Compute_jitter :
*	Convolve the PSF with a jitter function.
*
*  Inputs :
*	Crit_psf   :  Complex array containing critically sampled PSF.
*	dim        :  dimension of Crit_psf (dim x dim).
*	wavelength :  wavelength of PSF in microns.
*	jitter     :  RMS jitter in milliarceconds.
*			
*----------------------------------------------------------------------*/
void Compute_jitter( complex **crit_psf, int dim, float wavelength, 
		     float jitter_major, float jitter_minor, float jitter_angle )
{
	float	c, v, radius;
	int	x, y;


	radius = dim / 2;

	cos_angle = cos(angle * M_PI / 180);
	sin_angle = sin(angle * M_PI / 180);

	/* Multiply the OTF (FFT of PSF) by the jitter function */

	cx = M_PI * jitter_major * 4.848e-9 * DIAM_HST_MICRONS / wavelength;
	cy = M_PI * jitter_minor * 4.848e-9 * DIAM_HST_MICRONS / wavelength;

	for ( y = 0; y < dim; ++y )
	{
		if ( y < dim/2 )
			yj0 = y / radius;
		else
			yj0 = (dim - y) / radius;

		for ( x = 0; x < dim; ++x )
		{
			if ( x < dim/2 )
				xj0 = x / radius;
			else
				xj0 = (dim - x) / radius;

			xj0 = xj0 * cx*cx;
			yj0 = yj0 * cy*cy;
			xj = xj0*cos_angle - yj0*sin_angle;
			yj = xj0*sin_angle + yj0*cos_angle;

			v = sqrt(xj*xj + yj*yj);
			v = exp( -2.0 * v * v );
			crit_psf[y][x].r *= v;
			crit_psf[y][x].i *= 0.0;
		}
	}

} /* Compute_jitter */

