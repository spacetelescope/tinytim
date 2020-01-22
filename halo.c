/* File : halo.c
 *
 * Contents :
 *	Halo_f1042m : Add halo to WFPC2 F1042M PSFs
 *	Acs_pixel_area : Compute area of distorted ACS pixel
 *	Hrc_halo : Add ACS HRC red halo to PSF
 *
 * Written by John Krist, November 2000
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*----------------------------------------------------------------------------
*  Routine : Halo_f1042m 
*
*  Add WFPC2 F1042M filter halo to Image.  Image must be flux normalized
*  prior to calling this routine.  Subsamp is the ratio of the PSF pixel
*  size to the normal detector pixel size (ie. subsamp=0.5 if 2x smaller
*  than normal detector pixels).  Image is dim x dim pixels.
*
*  The halo equation is defined in terms of detector pixels (not arcsec).
*---------------------------------------------------------------------------*/
void Halo_f1042m( float **image, int dim, float subsamp )
{
	int	x, y, xc, yc;
	float	area, r, total, xr, yr;

	xc = dim / 2;
	yc = dim / 2;

	area = subsamp * subsamp;

	/* normalize total PSF flux to 0.805 */

	total = 0.0;
	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
			total += image[y][x];

	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
			image[y][x] *= 0.805 / total ;

	/* add in halo */

	for ( y = 0; y < dim; ++y )
	{
		yr = y - yc;
		for ( x = 0; x < dim; ++x )
		{
			xr = x - xc;
			r = sqrt(xr*xr + yr*yr) * subsamp;

			if ( r > 26.1 )
				image[y][x] += 0.3703*area * pow(r,-3.0);
			else if ( r > 5.7 )
				image[y][x] += 7.475e-4*area * pow(r,-1.1);
			else if ( r > 2.2 )
				image[y][x] += 0.0143*area * pow(r,-2.8);
			else
				image[y][x] += 0.0143*area * pow(2.2,-2.8);  
		}
	}

} /* Halo_f1042m */

/*---------------------------------------------------------------------------*/
float Acs_pixel_area( int x, int y, int camera )
{
	/* returns area of pixel at x,y in square arcsec */
 
	float	area, v2[4], v3[4];
	int	i, j;

	XY_to_V2V3( x-0.5, y-0.5, &v2[0], &v3[0], camera );
	XY_to_V2V3( x+0.5, y-0.5, &v2[1], &v3[1], camera );
	XY_to_V2V3( x+0.5, y+0.5, &v2[2], &v3[2], camera );
	XY_to_V2V3( x-0.5, y+0.5, &v2[3], &v3[3], camera );
	
	area = 0.0;

	for (i = 0; i < 4; ++i) 
   	{ 
      		j = (i==3) ? 0 : i+1;
      		area += (v3[j] + v3[i]) * (v2[j] - v2[i]) / 2; 
   	}

	return( area );
}

/*--------------------------------------------------------------------*/
void Hrc_halo( float **psf, float wavelength, int psf_x, int psf_y )
{
	float	a, b, f, r, t, xd, yd, v2c, v3c, v2, v3;
	int	x, y;

	v2v3( Pars.chip, psf_x, psf_y, &v2c, &v3c );  /* pupil v2,v3 */
	v2c = -v2c;	/* field coordinates are -pupil coordinates */
	v3c = -v3c;

	b = -2.36 * wavelength + 42.0; 

	a = Pars.psf_scale*Pars.psf_scale / Acs_pixel_area( psf_x, psf_y, ACS_HRC );
	f = a * 2.0 * pow(wavelength-0.45, 3.0) / 
		(-1181.62 * wavelength + 11036.7);

	for ( y = 0; y < Pars.int_dim; ++y )
	{
		v3 = (y - Pars.int_dim/2) * Pars.psf_scale + v3c;

		for ( x = 0; x < Pars.int_dim; ++x )
		{
			/* for HRC, +v2 is in -x direction */
			v2 = -(x - Pars.int_dim/2) * Pars.psf_scale + v2c;

			/* convert v2,v3 to distorted pixel coordinates */

			V2V3_to_XY( &xd, &yd, v2, v3, ACS_HRC );
			xd = xd - psf_x;
			yd = yd - psf_y;

			r = sqrt( xd*xd + yd*yd );
			t = -r / b;
			if ( t < -50.0 )
				t = -50.0;
			psf[y][x] = psf[y][x] * (1-f) + exp(t) * f;
		}
	}

} /* Hrc_halo */

