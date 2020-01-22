/* File     :  inter.c
 *
 * Contents :
 *	Interpolate  :  Bilinear interpolation (used by CreateMirrorMap)
 *
 * Author   :  John Krist (STScI)
 * Date     :  January 1992
 *
 * Modifications :
 *	J. Krist - April 1992
 *	  Modified Interpolate to use separate x and y starting values
 *	  for the images.
 */


#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*-------------------------------------------------------------------------
*  Interpolate :
*	Bilinear interpolation of a mirror map.
*
*  Inputs:
*	old_image    :  Image to be expanded/shrunk.
*	old_dim      :  old_image dimension (old_dim x old_dim).
*	old_start_x   :  Starting x of old_image in unit circle units from center.
*	old_start_y   :  Starting y of old_image in unit circle units from center.
*	old_spacing  :  Spacing of old_image in unit circle units.
*	new_image    :  Image containing interpolated image.
*	new_dim      :  Dimension of new_image.
*	new_start_x   :  Starting x of new_image, in unit circle units from center.
*	new_start_y   :  Starting y of new_image, in unit circle units from center.
*	new_spacing  :  Spacing of new_image in unit circle units.
*
*  Returns :
*	The new image is placed in new_image.
*--------------------------------------------------------------------------*/
void Interpolate( float **old_image, int old_dim, float old_start_x, float old_start_y, float old_spacing,
	          float **new_image, int new_dim, float new_start_x, float new_start_y, float new_spacing )
{
	float	diff_x[4000], old_pixel_x[4000], diff_y[4000], old_pixel_y[4000]; 
	float   old_end_x, old_end_y, t, u;
	int	xa, ya, x1, x2, y1, y2, x, y;


	old_end_x = old_start_x + old_spacing * (old_dim - 1);
	old_end_y = old_start_y + old_spacing * (old_dim - 1);

	x1 = 0;
	x2 = new_dim - 1;
	y1 = 0;
	y2 = new_dim - 1;

	/* Create lookup tables of differences and pixel values for a  *
  	 * given index.						       */

	for ( x = 0; x < new_dim; ++x )
	{
		t = x * new_spacing + new_start_x;
		if ( t < old_start_x )
			x1 = x + 1;
		else if ( t <= old_end_x )
		{
			x2 = x;
			u = (t - old_start_x) / old_spacing;
			old_pixel_x[x] = (int)u;
			diff_x[x] =  u - old_pixel_x[x];
		}

		t = x * new_spacing + new_start_y;
		if ( t < old_start_y )
			y1 = x + 1;
		else if ( t <= old_end_y )
		{
			y2 = x;
			u = (t - old_start_y) / old_spacing;
			old_pixel_y[x] = (int)u;
			diff_y[x] =  u - old_pixel_y[x];
		}
	}

	/* Do the bilinear interpolation */

	for ( y = y1; y <= y2; ++y )
	{
		ya = old_pixel_y[y];
		u = diff_y[y];
		for ( x = x1; x <= x2; ++x )
		{
			xa = old_pixel_x[x];
			t = diff_x[x];
			new_image[y][x] = (1. - t) * (1. - u) * old_image[ya][xa]
				+ t * (1. - u ) * old_image[ya][xa+1]
				+ t * u * old_image[ya+1][xa+1]
				+ (1. - t) * u * old_image[ya+1][xa];
		}
	}

}  /* Interpolate */	
