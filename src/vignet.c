/* File     :  vignet.c
 *
 * Contents :
 *	WfcVignette :  Multiply WFC (WFPC I) aperture by a circle with
 *		a blurred edge to simulate vignetting. Kludge.
 *
 * Author  : John Krist (STScI)
 * Date    : January 1993
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

#define Sqr(x)  ((x)*(x))

/*-------------------------------------------------------------------------*/
void Wfc_vignette( float **image, int nx, int ny )
{
	float	inner_radius, radius, xc, yc, *row, ysqr, c, d;
	int     x, y, x1, x2, y1, y2;

	radius = Pars.wfc_vig_radius * nx / 4;
	xc = Pars.wfc_vig_xc * nx / 4 + nx / 2;
	yc = Pars.wfc_vig_yc * nx / 4 + nx / 2;	
	inner_radius = Pars.wfc_vig_inside * radius;

	/* If inner_radius is 0.0, assume that you do not want to draw *
 	 * the vignetting.					      */

	if ( inner_radius <= 0.0 )
		return;

	/* The routine works like this:  If a pixel is outside of the *
	 * radius, it is set to 0.0; all pixels within some percent   *
	 * of the radius are kept as they are.  Pixels between the    *
	 * inner and outer radii are multiplied by a linear scale,    *
	 * going from 1.0 to 0.0 as the distance from the center      *
	 * increases.  In effect, this causes a "blurry" edge.        */

	/* Determine top and bottom rows which contain the circle */

	y1 = (int)(yc - radius);
	y2 = (int)(yc + radius) + 1;
	if ( y1 < 0 )
		y1 = 0;
	if ( y2 >= ny )
		y2 = ny - 1;

	/* Fill up rows before circle with zeroes */

	for ( y = 0; y < y1; ++y )
		for ( x = 0; x < nx; ++x )
			image[y][x] = 0.0;

	/* Fill up rows after circle with zeroes */

	for ( y = y2 + 1; y < ny; ++y )
		for ( x = 0; x < nx; ++x )
			image[y][x] = 0.0;

	for ( y = y1; y <= y2; ++y )
	{
		row = image[y];

		/* Determine the extent of the circle in this row */

		if ( y <= yc )
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y+1-yc) )) + 1;
		else
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y-1-yc) )) + 1;

		x1 = xc - x2;
		x2 = xc + x2;

		if ( x1 < 0 )
			x1 = 0;
		else if ( x2 >= nx )
			x2 = nx - 1; 

		/* Fill preceding and trailing columns with zeroes */

		for ( x = 0; x < x1; ++x )
			row[x] = 0.0;
		for ( x = x2 + 1; x < nx; ++x )
			row[x] = 0.0;

		ysqr = Sqr(y - yc);

		for ( x = x1; x <= x2; ++x )
		{
		     /* If pixel is greater than radius+2 from the center,  *
		      * the set it to 0.0; if it is less than radius-2 from *
		      * the center, then don't change it; if it is in       *
		      * between, subpixelate the pixel to determine how     *
		      * much of the circle goes through it, and multiply    *
		      * the pixel value appropriately.                      */

		     d = sqrt( Sqr(x-xc) + ysqr );

		     if ( d > radius )
			  row[x] = 0.0;
		     else if ( d >= inner_radius )
		     {
/*
			  c = 1.0 - (d - inner_radius) / (radius - inner_radius);
*/
			  c = (d - inner_radius) / (radius - inner_radius);
			  c = 1.0 - Sqr(c);
			  if ( c < 0.0 )
				c = 0.0;
			  else if ( c > 1.0 )
				c = 1.0;
			  row[x] *= c;
		     }
		}
	}

}  /* Wfc_vignette */
