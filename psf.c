/*  
 *  File     :   psf.c
 *
 *  Contents :   
 *	Compute_crit_psf  :  Compute critically sampled PSF
 *
 *  Author   :   John Krist (Space Telescope Science Institute)
 *
 *  Date     :   January 1992
 *
 *  Modifications :
 *	August 1994 - JEK
 *		Modified routines to use complex structures.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tinytim.h"


/*--------------------------------------------------------------------------
*  Compute_crit_psf :
*       Compute the critically sampled PSF from the pupil function (pupil)
*       and the optical path difference function (Opd).
*
*  Inputs :
*       pupil  :  Image array of dim by dim complex pixels containing the pupil
*                    function.
*       dim    :  dimension of pupil (dim x dim).
*
*  Outputs :
*       Replaces "pupil" with the critically sampled PSF.
*--------------------------------------------------------------------------*/
void Compute_crit_psf( complex **pupil, int dim )
{
	int     x, y;
	float   t, u, v, two_pi;

	two_pi = M_PI * 2.0;

	/* Compute  pupil.r * exp( 2*Pi*Func ) where Func = 0 + (pupil.i)i. *
	 * The complex exponential is replaced with its sin, cos terms.     */

	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
		{
			t = pupil[y][x].i * two_pi;
			pupil[y][x].i = sin(t) * pupil[y][x].r;
			pupil[y][x].r *= cos(t);
		}

	/* Compute the amplitude spread function (ASF) */

	fft2d( pupil[0], dim, 1 );

	/* Compute the PSF as the modulus squared of the ASF, place it  * 
	 * in the pupil array.                                          */

	for ( y = 0; y < dim; ++y )
	{
		for ( x = 0; x < dim; ++x )
		{
			t = pupil[y][x].r;
			u = pupil[y][x].i;
			v = sqrt( t*t + u*u );
			pupil[y][x].r = v * v;
			pupil[y][x].i = 0.0;
		}
	}

} /* Compute_crit_psf */

