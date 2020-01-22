/* File     :  pupil.c
 *
 * Contents :
 *      Draw_ellipse        : Draw an elliptical aperture (for STIS Lyot stop)
 *	Smooth	  	    : Smooth image using boxcar averaging
 *	Dark_inside_circle  : Draw a dark circle into the pupil function
 *				using subpixellation.
 *	Dark_outside_circle : Mask out everything the the pupil function
 *				outside of a circle.
 *	Clear_inside_circle : Draw a clear circle into the pupil function
 *				using subpixellation.
 *	Draw_rectangle      : Draw a subpixellated dark rectangle into the
 *				pupil function.
 *      Dark_rot_rect       : Draw a subpixellated rotated dark rectangle into the
 *              		pupil function using 
 *
 * Author  : John Krist (STScI)
 * Date    : January 1992
 *
 * Comments :  Okay, I know a lot of this code is redundant.  What do you
 *		want for free?
 *
 * Modifications:
 *	J. Krist - March 1992
 *	   Changed DarkInsideCircle and ClearInsideCircle to
 *	   do nothing if the circle radius is 0.0.
 *
 *      J. Krist - April 1992
 *	   Changed routines which draw circles to more accurately compute
 *	   the subpixel distances. 
 *	   Set the SUB_FACTOR to 10, which is what TIM uses.
 *
 *      R. N. Hook - October 1996
 *         Added DarkRotRect routine which makes dark rotated rectangles
 *         using the "Boxer" code. This was needed for the NICMOS pad
 *         cover obscurations.
 *
 *	J. Krist - April 1998
 *	   Added DrawEllipse and Smooth for STIS Lyot stop
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"
#include <stdlib.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include <memory.h>
#include <string.h>


#define SUB_FACTOR  8
#define Sqr(x)  ((x)*(x))

/*-------------------------------------------------------------------------*/
void Draw_ellipse( float **image, int nx, int ny, float xc, float yc, 
		   float x_radius, float y_radius )
{
	int     x, y;
	float	xx, yy;

	/* NOTE: ellipse is not antialiased, since it is assumed that it *
	 * will be smoothed (used for STIS Lyot stop)			 */

	for ( y = 0; y < ny; ++y )
	{
	   for ( x = 0; x < nx; ++x )
	   {
		xx = (x - xc) / x_radius;
		yy = (y - yc) / y_radius;
		if ( sqrt(xx*xx + yy*yy) > 1.0 )
			image[y][x] = 0.0;
		else
			image[y][x] = 1.0;
	   }
	}

}  /* Draw_ellipse */

/*----------------------------------------------------------------------------*/
void Smooth( float **image, int nx, int ny, int box_size )
{
	float	**buffer, *tempptr, Npixels;
	int	i, j, x, y;

	/* box_size/2 pixel perimeter around image is not modified */

	/* box_size must be odd-valued */

	Npixels = box_size * box_size;

	buffer = (float **)malloc( box_size * sizeof(float *) );

	for ( y = 0; y < box_size; ++y )
		buffer[y] = (float *)malloc( nx * sizeof(float) );

	for ( y = 1; y < box_size; ++y )
	{
		for ( x = 0; x < nx; ++x )
			buffer[y][x] = image[y-1][x];
	}

	for ( y = box_size/2; y < ny-box_size/2; ++y )
	{
	   tempptr = buffer[0];
	   for ( j = 0; j < box_size-1; ++j )
		buffer[j] = buffer[j+1];
	   buffer[box_size-1] = tempptr;
	   for ( i = 0; i < nx; ++i )
		buffer[box_size-1][i] = image[y+box_size/2][i];

	   for ( x = box_size/2; x < nx-box_size/2; ++x )
	   {
		for ( j = 0; j < box_size; ++j )
		   for ( i = x-box_size/2; i <= x+box_size/2; ++i )
			image[y][x] += buffer[j][i];
		image[y][x] /= Npixels;
	   }
	}
		
	for ( j = 0; j < box_size; ++j )
		free( buffer[j] );
	free( buffer );

} /* Smooth */

/*----------------------------------------------------------------------------*/
void Dark_inside_circle( float **image, int nx, int ny, float xc, float yc, float radius )
{
	float   *row, radius_min_2, value_table[200];
	float   xsub, ysub, xsub2, ysub2, sub_inc, ysqr;
	int     x, y, x1, x2, y1, y2, sub_pixels_set;


	if ( radius < 1.e-6 )
		return;

	radius_min_2 = radius - 2.0;
	sub_inc = 1.0 / SUB_FACTOR;

	/* Fill lookup table which converts number of subpixels set to  *
	 * percentage of the pixel set.                                 */

	for ( x = 0; x <= SUB_FACTOR*SUB_FACTOR; ++x )
		value_table[x] = x * (sub_inc*sub_inc);

	y1 = (int)(yc - radius); 
	y2 = (int)(yc + radius) + 1;
	if ( y1 < 0 )
		y1 = 0;
	if ( y2 >= ny )
		y2 = ny - 1;

	for ( y = y1; y <= y2; ++y )
	{
		row = image[y];
		if ( y <= yc )
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y+1-yc) )) + 1;
		else
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y-1-yc) )) + 1;

		x1 = xc - x2;
		x2 = x2 + xc;

		if ( x1 < 0 )
			x1 = 0;
		else if ( x2 >= nx )
			x2 = nx - 1; 

		ysqr = Sqr(y - yc);

		for ( x = x1; x <= x2; ++x )
		{
		     if ( sqrt( Sqr(x-xc) + ysqr ) > radius_min_2 )
		     {
			  sub_pixels_set = 0;
			  for ( ysub = y; ysub < y+0.999; ysub += sub_inc )
			  {
				if ( ysub < yc )
					ysub2 = Sqr(ysub + sub_inc - yc);
				else
					ysub2 = Sqr(ysub - yc);

				for ( xsub = x; xsub < x+0.999; xsub+= sub_inc )
				{
					if ( xsub < xc )
						xsub2 = Sqr(xsub + sub_inc - xc);
					else
						xsub2 = Sqr(xsub - xc);

					if ( sqrt(xsub2 + ysub2) > radius )
						++sub_pixels_set;
				}
			  }
			  row[x] *= value_table[sub_pixels_set];
		     }
		     else
			  row[x] = 0.0;
		}
	}

}  /* Dark_inside_circle */

/*-------------------------------------------------------------------------*/
void Dark_outside_circle( float **image, int nx, int ny, float xc, float yc, float radius )
{
	float   xsub, ysub, xsub2, ysub2, sub_inc;
	float	*row, value_table[200];
	float   radius_min_2, radiusPlus2, ysqr, d;
	int     x, y, x1, x2, y1, y2, sub_pixels_set;


	sub_inc = 1.0 / SUB_FACTOR;   /* Subpixel size along one axis */

	/* Fill lookup table which converts number of subpixels set to  *
	 * percentage of the pixel set.                                 */

	for ( x = 0; x <= SUB_FACTOR*SUB_FACTOR; ++x )
		value_table[x] = x * (sub_inc*sub_inc);

	radius_min_2 = radius - 2.0;
	radiusPlus2 = radius + 2.0;

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

		     if ( d > radiusPlus2 )
			  row[x] = 0.0;
		     else if ( d >= radius_min_2 )
		     {
			  sub_pixels_set = 0;
			  for (ysub=y; ysub < y + 0.999;ysub += sub_inc )
			  {
				if ( ysub < yc )
					ysub2 = Sqr(ysub + sub_inc - yc);
				else
					ysub2 = Sqr(ysub - yc);

				for (xsub=x;xsub<x+0.999;xsub += sub_inc)
				{
					if ( xsub < xc )
						xsub2 = Sqr(xsub + sub_inc - xc);
					else
						xsub2 = Sqr(xsub - xc);

					if ( sqrt(xsub2 + ysub2) <= radius )
						++sub_pixels_set;
				}
			  }
			  row[x] *= value_table[sub_pixels_set];
		     }
		}
	}

}  /* Dark_outside_circle */

/*-------------------------------------------------------------------------*/
void Clear_inside_circle( float **image, int nx, int ny, float xc, float yc, float radius )
{
	float   xsub, ysub, xsub2, ysub2, sub_inc;
	float   *row, radius_min_2, ysqr, value_table[200];
	int     x, y, x1, x2, y1, y2, sub_pixels_set;


	if ( radius < 1.e-6 )
		return;

	radius_min_2 = radius - 2.0;
	sub_inc = 1.0 / SUB_FACTOR;

	/* Fill lookup table which converts number of subpixels set to  *
	 * percentage of the pixel set.                                 */

	for ( x = 0; x <= SUB_FACTOR*SUB_FACTOR; ++x )
		value_table[x] = x * (sub_inc*sub_inc);

	y1 = (int)(yc - radius); 
	y2 = (int)(yc + radius) + 1;
	if ( y1 < 0 )
		y1 = 0;
	if ( y2 >= ny )
		y2 = ny - 1;

	for ( y = y1; y <= y2; ++y )
	{
		row = image[y];
		if ( y <= yc )
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y+1-yc) )) + 1;
		else
			x2 = (int)(sqrt( Sqr(radius) - Sqr(y-1-yc) )) + 1;

		x1 = xc - x2;
		x2 = x2 + xc;

		if ( x1 < 0 )
			x1 = 0;
		else if ( x2 >= nx )
			x2 = nx - 1; 

		ysqr = Sqr(y - yc);

		for ( x = x1; x <= x2; ++x )
		{
		     if ( sqrt( Sqr(x-xc) + ysqr ) > radius_min_2 )
		     {
			  sub_pixels_set = 0;
			  for (ysub=y; ysub < y+0.999; ysub += sub_inc )
			  {
				if ( ysub < yc )
					ysub2 = Sqr(ysub + sub_inc - yc);
				else
					ysub2 = Sqr(ysub - yc);

				for (xsub=x; xsub<x+0.999;xsub+=sub_inc )
				{
					if ( xsub < xc )
						xsub2 = Sqr(xsub + sub_inc - xc);
					else
						xsub2 = Sqr(xsub - xc);

					if ( sqrt(xsub2 + ysub2) <= radius )
						++sub_pixels_set;
				}
			  }
			  row[x] = value_table[sub_pixels_set];
		     }
		     else
			  row[x] = 1.0;
		}
	}

}  /* Clear_inside_circle */

/*-------------------------------------------------------------------------*/
void Dark_rectangle( float **image, int nx, int ny, float xc, float yc, 
		     float x_length, float y_length )
{
	float   left, right, top, bottom;
	int     i, j, left_pixel, right_pixel, top_pixel, bottom_pixel;


	/* Determine how much of the edge pixels are filled */

	if ( x_length < 1.e-5 || y_length < 1.e-5 )
		return;

	left = xc - x_length / 2.0;  /* distance of left edge from center */
	left_pixel = (int)left;
	if ( left > 0.0 )
		left = left - left_pixel;
	else
		left = 0.0;

	right = xc + x_length / 2.0;
	right_pixel = (int)right;
	if ( right < nx )
		right = 1.0 - right + right_pixel;
	else
		right = 0.0;

	bottom = yc - y_length / 2.0;
	bottom_pixel = (int)bottom;
	if ( bottom > 0.0 )
		bottom = bottom - bottom_pixel;    
	else
		bottom = 0.0;

	top = yc + y_length / 2.0; 
	top_pixel = (int)top;
	if ( top < ny )
		top = 1.0 - top + top_pixel;
	else 
		top = 0.0;

	for ( i = left_pixel; i <= right_pixel; ++i )
	{
		image[bottom_pixel][i] *= bottom;
		image[top_pixel][i] *= top;
	}

	for ( j = bottom_pixel; j <= top_pixel; ++j )
	{
		image[j][left_pixel] *= left;
		image[j][right_pixel] *= right;
	}

	for ( j = bottom_pixel + 1; j < top_pixel; ++j )
		for ( i = left_pixel + 1; i < right_pixel; ++i )
			image[j][i] = 0.0;

}  /* Dark_rectangle */

/*----------------------------------------------------------------------*/
float Sgarea( float x1, float y1, float x2, float y2 )
{
	/* To calculate area under a line segment within unit  *
	 * square at origin. This is used by BOXER             */

      	float m, c, dx;
      	float xlo, xhi, ylo, yhi, xtop;

      	dx = x2 - x1;

	/* Trap vertical line */

      	if (dx == 0.0) 
		return (0.0);

	/* Order the two input points in x */

      	if (x1 < x2)  
	{
         	xlo = x1;
         	xhi = x2;
        }
 	else 
	{
         	xlo = x2;
         	xhi = x1;
        }

	/* And determine the bounds ignoring y for now */

      	if (xlo >= 1.0 || xhi <= 0.0) 
		return (0.0); 

      	if ( xlo < 0.0) 
		xlo=0.0;
      	if ( xhi > 1.0) 
		xhi=1.0;

	/* Now look at y basic info about the line y = mx + c */

      	m = (y2 - y1) / dx;
      	c = y1 - m * x1;
      	ylo = m * xlo + c;
      	yhi = m * xhi + c;

	/* Trap segment entirely below axis */

      	if (ylo <= 0.0 && yhi <= 0.0) 
		return(0.0);

	/* Adjust bounds if segment crosses axis (to exclude anything below axis) */

      	if (ylo < 0.0)  
	{
         	ylo = 0.0;
         	xlo = -c/m;
      	}

      	if (yhi < 0.0)  
	{
         	yhi = 0.0;
         	xhi = -c/m;
      	}

	/* There are four possibilities: both y below 1, both y above 1
   	 * and one of each. */

      	if (ylo >= 1.0 && yhi >= 1.0)  
        {
		/* Line segment is entirely above square */

         	if (dx < 0.0) 
            		return(xlo - xhi);
         	else
            		return(xhi - xlo);
        }

      	if (ylo <= 1.0 && yhi <= 1.0)  
        {
		/* Segment is entirely within square */

         	if (dx < 0.0) 
            		return(0.5 * (xlo-xhi) * (yhi+ylo));
         	else
            	return(0.5 * (xhi-xlo) * (yhi+ylo));
        }

	/* otherwise it must cross the top of the square */

      	xtop = (1.0 - c) / m;

	/* Check for a problem */

      	if (xtop < xlo || xtop > xhi) 
       		printf("! Warning, box overlap calculation problem\n"); 

      	if (ylo < 1.0)  
	{
         	if (dx < 0.0)  
            		return(-(0.5 * (xtop-xlo) * (1.0+ylo) + xhi - xtop));
         	else
            		return(0.5 * (xtop-xlo) * (1.0+ylo) + xhi - xtop);
      	}

      	if (dx < 0.0) 
         	return(-(0.5 * (xhi-xtop) * (1.0+yhi) + xtop-xlo));
      	else
         	return(0.5 * (xhi-xtop) * (1.0+yhi) + xtop-xlo);
}

/*-------------------------------------------------------------------------
* Boxer :
*	Compute area of box overlap
*
* Calculate the area common to input clockwise polygon x(n), y(n) with
* square (is, js) to (is+1, js+1).
* This version is for a quadrilateral.
*
* W.B. Sparks STScI 2-June-1990.
* Phil Hodge        20-Nov-1990  Change calling sequence; single precision.
* Richard Hook ECF  24-Apr-1996  Change coordinate origin
*                                so that centre of pixel has integer position
* Richard Hook ECF  8-Oct-1996   Convert from F77 to C
*------------------------------------------------------------------------*/
float Boxer ( int is, int js, float *x, float *y )
{
      float px[4], py[4], sum;
      int   i;

	/* Set up coords relative to unit square at origin            *
	 * Note that the +0.5s were added when this code was modified */

	for ( i = 0; i < 4; ++i ) 
	{
         	px[i] = x[i] - (float)is + 0.5;
         	py[i] = y[i] - (float)js + 0.5;
        } 
 
	/* For each line in the polygon (or at this stage, input quadrilateral)
   	 * calculate the area common to the unit square (allow negative area for
   	 * subsequent `vector' addition of subareas). */

      	sum = 0.0;
      	for ( i = 0; i < 3; ++i ) 
         	sum +=  Sgarea(px[i], py[i], px[i+1], py[i+1]);

      	sum += Sgarea(px[3], py[3], px[0], py[0]);

      	return(sum);
}

/*-------------------------------------------------------------------------*/
void Dark_rot_rect( float **image, int nx, int ny, float xc, float yc, 
		    float x_length, float y_length, float theta )
{
	float   x[4], y[4], phi, r;
        float   xmin, xmax, ymin, ymax;
        int     ist, iend, jst, jend;
	int     i, j;

        /* Work out where the four corners lie */

	if ( y_length != 0.0 || x_length != 0.0 )
        	phi = atan2( y_length, x_length );
	else
		phi = 0.0;

        r = 0.5 * sqrt( x_length*x_length + y_length*y_length);

        x[0] = xc + r * cos( phi - theta);
        x[3] = xc + r * cos( M_PI -phi - theta);
        x[2] = xc + r * cos( phi - theta - M_PI);
        x[1] = xc + r * cos( -phi - theta);
        y[0] = yc + r * sin( phi - theta);
        y[3] = yc + r * sin( M_PI -phi - theta);
        y[2] = yc + r * sin( phi - theta - M_PI);
        y[1] = yc + r * sin( -phi - theta);

        /* Work out the range of X and Y covered and clip
         * at the edges. */

        xmin = x[0];
        xmax = x[0];
        ymin = y[0];
        ymax = y[0];

        for ( i = 1; i < 4; ++i ) 
	{
            	if ( x[i] < xmin ) 
			xmin = x[i];
            	if ( x[i] > xmax ) 
			xmax = x[i];
            	if ( y[i] < ymin ) 
			ymin = y[i];
            	if ( y[i] > ymax ) 
			ymax = y[i];
        }

        ist = (int)(xmin - 1.0);
        iend = (int)(xmax + 1.0);
        jst = (int)(ymin - 1.0);
        jend = (int)(ymax + 1.0);

        if ( ist < 0 ) 
		ist = 0;
        if ( iend > nx-1 ) 
		iend = nx - 1;
        if ( jst < 0 ) 
		jst = 0;
        if ( jend> ny-1 ) 
		jend = ny - 1;

        /* Now loop over the region and work out the overlaps */

        for ( j = jst; j < jend; ++j )
            	for ( i = ist; i < iend; ++i )
                	image[j][i] *= (1.0 - Boxer( i, j, x, y));

}  /* Dark_rot_rect */

