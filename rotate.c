/*  File     :  rotate.c
 *
 *  Contents :
 *      Rotateimage  :  Rotate an image
 *  	Rotate : Rotate, scale, and shift an image
 *
 *  Author   :  John Krist (STScI)
 *  Date     :  January 1992
 *
 *  Modifications :
 *	JEK - March 1994
 *		Added Rotate routine.
 *
 *	JEK - April 1994
 *		Fixed bug in Rotate routine (allocating memory without
 *		zeroing it).
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

#define  MAX_ARRAY   4096


/*----------------------------------------------------------------
*  Rotate_image  :  Rotate a 2-d image through an arbitrary
*                   angle.  Uses bilinear interpolation.
*
*       image   =  2-d array of the form image[ny][nx].
*       nx, ny  =  image dimensions.
*       xc, yc  =  Center of rotation.
*       theta   =  Clockwise rotation angle (in degrees, and
*                  assuming pixel (0,0) is in the lower left).
*
*       Rotateimage returns a pointer to a newly created array
*       which contains the rotated image.
*----------------------------------------------------------------*/
float **Rotate_image( float **image, int nx, int ny, float xc, float yc, float theta )
{
	float   **temp;
	float   x_cos_theta[MAX_ARRAY], x_sin_theta[MAX_ARRAY];
	float   y_cos_theta[MAX_ARRAY], y_sin_theta[MAX_ARRAY];
	float   old_x, old_y, t, u, z1, z2, z3, z4;
	int     x, y, x1, x2, y1, y2;


	temp = Alloc_image( nx, ny );

	theta = theta * M_PI / 180.0;   /* Convert theta to radians */

	/* Load transform constants into lookup tables */

	for ( x = 0; x < nx; ++x )
	{
		x_cos_theta[x] = (x - xc) * cos(theta);
		x_sin_theta[x] = (x - xc) * sin(theta);
	}

	for ( y = 0; y < ny; ++y )
	{
		y_cos_theta[y] = (y - yc) * cos(theta);
		y_sin_theta[y] = (y - yc) * sin(theta);
	}

	for ( y = 0; y < ny; ++y )
	{
		for ( x = 0; x < nx; ++x )
		{
			old_x = x_cos_theta[x] - y_sin_theta[y] + xc;
			if ( old_x < 0.0 || old_x >= nx-1 )
				continue;
			old_y = x_sin_theta[x] + y_cos_theta[y] + yc;
			if ( old_y < 0.0 || old_y >= ny-1 )
				continue;
	
			x1 = (int)old_x;
			x2 = x1 + 1;
			y1 = (int)old_y;
			y2 = y1 + 1;

			z1 = image[y1][x1];
			z2 = image[y1][x2];
			z3 = image[y2][x2];
			z4 = image[y2][x1];

			t = old_x - x1;
			u = old_y - y1;

			temp[y][x] = (1.0-t)*(1.0-u)*z1 + t*(1.0-u)*z2 +
					  t*u*z3 + (1.0-t)*u*z4;
		}
	}

	return( temp );

}  /* Rotate_image */

/*----------------------------------------------------------------
*  Rotate   :  Rotate a 2-d image through an arbitrary angle with scaling and
*               shifting.  Uses bilinear interpolation.  The process proceeds
*               effectively as follows :
*               1) Rotate image around xc, yc by theta degrees clockwise.
*               2) magnify the image by mag times.
*               3) Shift the image so that the center of rotation is located
*                    at new_xc, new_yc.
*
*       image   =  2-d array of the form image[ny][nx].
*       old_nx, old_ny  =  Input image dimensions.
*       xc, yc  =  Center of rotation in old image space.
*       new_nx, new_ny  =  Output image dimensions.
*       new_xc,
*       new_yc   =  Center of rotation in new image.
*       theta   =  Clockwise rotation angle (in degrees, and
*                  assuming pixel (0,0) is in the lower left).
*       mag     =  magnification
*
*       Rotate returns a pointer to a newly created array
*       which contains the new image.
*----------------------------------------------------------------*/
float **Rotate( float **image, int old_nx, int old_ny, float xc, float yc, 
		int new_nx, int new_ny, float new_xc, float new_yc, 
		float theta, float mag )
{
	float   **temp;
	float   x_cos_theta[MAX_ARRAY], x_sin_theta[MAX_ARRAY];
	float   y_cos_theta[MAX_ARRAY], y_sin_theta[MAX_ARRAY];
	float   old_x, old_y, t, u, xx, yy, z1, z2, z3, z4;
	int     x, y, x1, x2, y1, y2;


	temp = Alloc_image( new_nx, new_ny );

	theta = theta * M_PI / 180.0;   /* Convert theta to radians */

	/* Load transform constants into lookup tables */

	for ( x = 0; x < new_nx; ++x )
	{
		xx = (x - new_xc) / mag;
		x_cos_theta[x] = xx * cos(theta);
		x_sin_theta[x] = xx * sin(theta);
	}

	for ( y = 0; y < new_ny; ++y )
	{
		yy = (y - new_yc) / mag;
		y_cos_theta[y] = yy * cos(theta);
		y_sin_theta[y] = yy * sin(theta);
	}

	for ( y = 0; y < new_ny; ++y )
	{
		for ( x = 0; x < new_nx; ++x )
		{
			old_x = x_cos_theta[x] - y_sin_theta[y] + xc;
			if ( old_x < 0.0 || old_x >= old_nx-1 )
				continue;
			old_y = x_sin_theta[x] + y_cos_theta[y] + yc;
			if ( old_y < 0.0 || old_y >= old_ny-1 )
				continue;
	
			x1 = (int)old_x;
			x2 = x1 + 1;
			y1 = (int)old_y;
			y2 = y1 + 1;

			z1 = image[y1][x1];
			z2 = image[y1][x2];
			z3 = image[y2][x2];
			z4 = image[y2][x1];
	
			t = old_x - x1;
			u = old_y - y1;

			temp[y][x] = (1.0-t)*(1.0-u)*z1 + t*(1.0-u)*z2 +
					  t*u*z3 + (1.0-t)*u*z4;
		}
	}

	return( temp );

}  /* Rotate */

