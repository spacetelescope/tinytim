/*  File     :  cpupil.c
 *
 *  Contents :
 *    	Get_wfpc_shift    :  Compute WFPC obscuration shifts
 *      Draw_wfpc_spiders :  Draw WFPC spiders in pupil function
 *	Compute_pupil    :  Generate pupil function
 *
 *  Author   :  John Krist, STScI
 *  Date     :  January 1992
 *
 *  Modifications :
 *	J. Krist - April 1992
 *		Added code to handle WFPC II.
 *
 *	J. Krist - Aug 1992
 *		Removed DrawPC2Spiders.
 *
 *      R. Hook - Apr 1996
 *              Added code for NICMOS cold mask obscuration
 *
 *      J. Krist - June 1997
 *		Modified NICMOS
 *
 *	J. Krist - Sept 1997
 *		Added NICMOS cold mask offset parameters
 *
 *	J. Krist - Aug 1999
 *		Added NICMOS & OTA obscuration rotations
 *
 *      R. Hook & F. Stoehr - April 2008
 *              Added WFC3 IR channel support
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*----------------------------------------------------------------------------
*  Get_wfpc_shift :
*	Compute the center of the WFPC obscurations in (V3,V2) coordinates
*	for a given chip and detector position (from the TIM Users Manual).
*
*  Inputs :
*	chip   =  chip on which PSF is "imaged" : (1-8) or WFPC2_WFC_2-4 or
*			WFPC2_PC.
*	psf_x, psf_y =  Pixel coordinates of PSF, in detector coordinates.
*
*  Outputs:
*	wfpc_v3, wfpc_v2 :  V2,V3 coordinates, in unit circle units, which defines
*			 the center of the WFPC obscurations.
*----------------------------------------------------------------------------*/
static void Get_wfpc_shift( int chip, int psf_x, int psf_y, float *wfpc_v3, float *wfpc_v2 )
{
	float   i0, j0, xoff, yoff, xc, yc, xp, xm;


	i0 = Pars.obsc_center_x;
	j0 = Pars.obsc_center_y;

	/* Wide Field Camera and WFPC II WFC */

	if ( chip <= 4 || (chip >= WFPC2_WFC_2 && chip <= WFPC2_WFC_4))
	{
		xoff = 0.4000455 / 800.0; 
		yoff = 0.3880455 / 800.0;
		xc = xoff * ( psf_x - i0 - psf_y + j0 );
		yc = yoff * ( psf_x - i0 + psf_y - j0 );
	} 
	else if ( chip == WFPC2_PC ) 		/* WFPC2 Planetary Camera */
	{
		xc = 5.088e-4 * ( psf_x - i0 );
		yc = 5.088e-4 * ( psf_y - j0 );
	}
	else    /* WFPC1 Planetary Camera */
	{
		xoff = 0.1996364;
		yoff = 0.1936364;
		xp = xoff + yoff;
		xm = xoff - yoff;
		xc = (xp * (psf_x - i0) - xm*(psf_y - j0)) / (800.0 * sqrt(2.));
		yc = (-xm * (psf_x - i0) + xp*(psf_y - j0)) / (800.0 * sqrt(2.));
	}

	switch (chip)
	{
		case WFPC2_PC :
			  xc = -xc;
			  *wfpc_v2 = xc * cos(-M_PI/4.) - yc * sin(-M_PI/4.);
			  *wfpc_v3 = xc * sin(-M_PI/4.) + yc * cos(-M_PI/4.);
			  break;

		case 1 :  
		case 8 :
			  *wfpc_v2 = -xc;
			  *wfpc_v3 = yc;
			  break;

		case 2 :  
		case 5 :
		case WFPC2_WFC_2 :
			  *wfpc_v2 = yc;
			  *wfpc_v3 = xc;
			  break;

		case 3 :  
		case 6 :
		case WFPC2_WFC_3 :
			  *wfpc_v2 = xc;
			  *wfpc_v3 = -yc;
			  break;

		case 4 :  
		case 7 :
		case WFPC2_WFC_4 :
			  *wfpc_v2 = -yc;
			  *wfpc_v3 = -xc;
			  break;
	}

} /* Get_wfpc_shift */

/*-----------------------------------------------------------------------------
*  Draw_wfpc_spiders :
*	Draw the spiders from the WFPC repeater cameras into the pupil function.
*
*  Inputs:
*	pupil  :  Alloc_image-type array which contains the pupil function and
*		    into which the spiders will be drawn.
*	   dim :  dimensions of pupil (dim by dim)
*	wfpc_v2, wfpc_v3 :  Coordinates of WFPC obscuration centers in unit circle units.
*	chip   :  chip being used (1-8, WFPC2_WFC_2 to WFPC2_PC).
*----------------------------------------------------------------------------*/
static void Draw_wfpc_spiders( float **pupil, int dim, float wfpc_v2, float wfpc_v3, int chip )
{
	float   radius, width, length, center;


	radius = dim / 4.0;
	center = (int)(dim / 2) + 0.5;
	width = Pars.camera_spider_width * radius;

	/* Draw long and short spiders, in that order */

	switch (chip)
	{
		case 1: 
		case 7:
			Dark_rectangle( pupil, dim, dim, 
				radius * wfpc_v3 + center, 
				center, width, 2.2 * radius );
			length = (1.0 - wfpc_v3 + 0.1) * radius;
			Dark_rectangle( pupil, dim, dim, 
				radius*(wfpc_v3 + 1.)/2. + center, 
				radius * wfpc_v2 + center, 
				length, width );
			break;

		case WFPC2_PC :  /* short spider on other side */
			Dark_rectangle( pupil, dim, dim, 
				radius * wfpc_v3 + center, 
				center, width, 2.2 * radius );
			length = (1.0 + wfpc_v3 + 0.1) * radius;
			Dark_rectangle( pupil, dim, dim, 
				radius*(wfpc_v3 - 1.)/2. + center, 
				radius * wfpc_v2 + center, 
				length, width );
			break;

		case 2: 
		case 8:
		case WFPC2_WFC_2 :
			Dark_rectangle( pupil, dim, dim, center, 
				radius * wfpc_v2 + center, 
				2.2 * radius, width );
			length = (1.0 - wfpc_v2 + 0.1) * radius;
			Dark_rectangle( pupil, dim, dim, 
				radius * wfpc_v3 + center,
				radius * (wfpc_v2 + 1.) / 2. + center, 
				width, length );
			break;

		case 3: 
		case 5:
		case WFPC2_WFC_3 :
			Dark_rectangle( pupil, dim, dim, 
				radius * wfpc_v3 + center, 
				center, width, 2.2 * radius );
			length = radius * (wfpc_v3 + 1.0 + 0.1);
			Dark_rectangle( pupil, dim, dim, 
				radius * (wfpc_v3 - 1.) / 2. + center, 
				radius * wfpc_v2 + center,
				length, width );
			break;
		
		case 4: 
		case 6:
		case WFPC2_WFC_4 :
			Dark_rectangle( pupil, dim, dim, center, 
				radius * wfpc_v2 + center,
				2.2 * radius, width );
			length = radius * (wfpc_v2 + 1.0);
			Dark_rectangle( pupil, dim, dim,
				radius * wfpc_v3 + center,
				radius * (wfpc_v2 - 1.0) / 2. + center, 
				width, length );
			break;
	}

} /* Draw_wfpc_spiders */

/*----------------------------------------------------------------------------
*  Compute_pupil :
*	Create the pupil function array.
*
*  Inputs :
*       dim : dimension of pupil array (dim by dim)
*	psf_x, psf_y :  X,Y position on detector (ignored if no field dependence) 
*	pupil :  Pointer to image allocated with Alloc_image.
*
*  Returns:
*	Returns the pupil function in pupil.
*----------------------------------------------------------------------------*/
void Compute_pupil( int dim, int psf_x, int psf_y, float **pupil )
{
	int     i, x, y;
	float   radius, **newimage, wfpc_v2=-1e20, wfpc_v3=-1e20, center;
	float	nicmos_x_offset, nicmos_y_offset, **temp, **temp1;
	float	wfc3_x_offset, wfc3_y_offset;
        float   stis_x_offset, stis_y_offset;
	float	acs_x_offset, acs_y_offset;

	/* The pupil array starts off with pixel (0,0) in the lower left   *
	 * corner, with the V2,V3 origin at (dim/2+0.5,dim/2+0.5).  The    *
	 * V3 axis is horizontal, increasing to the right.  V2 is vertical *
	 * and increases upwards.   The pupil is later rotated to the      *
 	 * correct orientation for the specified camera.		   */

	center = dim / 2 + 0.5;
	radius = dim / 4.0;

	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
			pupil[y][x] = 1.0;

	/* Compute WFPC obscuration shifts */

	if ( (Pars.chip >= 1 && Pars.chip <= 8) ||
	     (Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4) )
		Get_wfpc_shift( Pars.chip, psf_x, psf_y, &wfpc_v3, &wfpc_v2 );

	/* Draw OTA aperture */

	Dark_outside_circle( pupil, dim, dim, center, center, radius );

	/* Draw OTA secondary */

	Dark_inside_circle( pupil, dim, dim, center, center, Pars.ota_secondary_radius * radius );

	/* Draw OTA spiders */

	Dark_rectangle( pupil, dim, dim, center, center, radius * 2. + 4.,
		Pars.ota_spider_width * radius );
	Dark_rectangle( pupil, dim, dim, center, center,
		Pars.ota_spider_width * radius, radius * 2. + 4. );

	/* Draw mirror pads */

	for ( i = 0; i < 3; ++i )
		Dark_inside_circle( pupil, dim, dim,
			Pars.ota_pad_v3[i] * radius + center,
			Pars.ota_pad_v2[i] * radius + center,
			Pars.ota_pad_radius[i] * radius );

	/* Draw WFPC obscurations */

	if ( (Pars.chip >= 1 && Pars.chip <= 8) ||
	     (Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4) )
	{
		/* Draw WFPC secondary obscuration */

		Dark_inside_circle( pupil, dim, dim,
				wfpc_v3 * radius + center,
				wfpc_v2 * radius + center,
				Pars.camera_secondary_radius * radius );

		/* Draw WFPC spiders, one full spider, one half spider */

	    	Draw_wfpc_spiders( pupil, dim, wfpc_v2, wfpc_v3, Pars.chip );

		/* If WFPC1 WFC (but not WFPC2 WFC), draw entrance aperture */

		if ( Pars.chip >= 1 && Pars.chip <= 4 )
		   Dark_outside_circle( pupil, dim, dim, wfpc_v3 * radius + center,
		       wfpc_v2 * radius + center, Pars.camera_entrance_radius * radius );
	}

	/* ACS coronagraph off-spot pupil */

	if ( Pars.chip == ACS_HRC_OFFSPOT )
	{
		acs_x_offset = Pars.acs_stop_shift * cos(Pars.acs_stop_angle*M_PI/180);
		acs_y_offset = Pars.acs_stop_shift * sin(Pars.acs_stop_angle*M_PI/180);

         	Dark_outside_circle( pupil, dim, dim, acs_x_offset*radius+center, acs_y_offset*radius+center, 
			Pars.acs_outer_radius * radius );
         	Dark_inside_circle( pupil, dim, dim, acs_x_offset*radius+center, acs_y_offset*radius+center, 
			Pars.acs_secondary_radius*radius );
         	Dark_outside_circle( pupil, dim, dim, -acs_x_offset*radius+center, -acs_y_offset*radius+center, 
			Pars.acs_outer_radius * radius );
         	Dark_inside_circle( pupil, dim, dim, -acs_x_offset*radius+center, -acs_y_offset*radius+center, 
			Pars.acs_secondary_radius*radius );

		for ( i = 0; i < 3; ++i )
			Dark_inside_circle( pupil, dim, dim,
				(Pars.ota_pad_v3[i]+acs_x_offset)*radius + center,
				(Pars.ota_pad_v2[i]+acs_y_offset)*radius + center,
				Pars.acs_pad_radius * radius );

		for ( i = 0; i < 3; ++i )
			Dark_inside_circle( pupil, dim, dim,
				(Pars.ota_pad_v3[i]-acs_x_offset)*radius + center,
				(Pars.ota_pad_v2[i]-acs_y_offset)*radius + center,
				Pars.acs_pad_radius * radius );

		Dark_rectangle( pupil, dim, dim, acs_x_offset*radius+center, acs_y_offset*radius+center, 
			radius * 2. + 4., Pars.acs_spider_width * radius );
		Dark_rectangle( pupil, dim, dim, acs_x_offset*radius+center, acs_y_offset*radius+center, 
			Pars.acs_spider_width * radius, radius * 2. + 4. );

		Dark_rectangle( pupil, dim, dim, -acs_x_offset*radius+center, -acs_y_offset*radius+center, 
			radius * 2. + 4., Pars.acs_spider_width * radius );
		Dark_rectangle( pupil, dim, dim, -acs_x_offset*radius+center, -acs_y_offset*radius+center, 
			Pars.acs_spider_width * radius, radius * 2. + 4. );
	}

        /* NICMOS cold mask obscuration */

	if ( (Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3)
	   || (Pars.chip >= NEWNICMOS_1 && Pars.chip <= NEWNICMOS_3) )
        {
		temp = Alloc_image( dim, dim );
		for ( y = 0; y < dim; ++y )
			for ( x = 0; x < dim; ++x )
				temp[y][x] = 1.0;

		nicmos_x_offset = Pars.nicmos_mask_x_off * radius + center;
		nicmos_y_offset = Pars.nicmos_mask_y_off * radius + center;

		/* Define inner, outer aperture radii */

         	Dark_outside_circle( temp, dim, dim, nicmos_x_offset, nicmos_y_offset, 
			Pars.camera_entrance_radius * radius );
         	Dark_inside_circle( temp, dim, dim, nicmos_x_offset, nicmos_y_offset, 
			Pars.camera_secondary_radius * radius );

		/* Draw NICMOS spiders (spiders are slighty non-perpendicular *
		 * from each other.					      */

		temp1 = Alloc_image( dim, dim );
		for ( y = 0; y < dim; ++y )
			for ( x = 0; x < dim; ++x )
				temp1[y][x] = 1.0;
         	Dark_rectangle( temp1, dim, dim, nicmos_x_offset, nicmos_y_offset, 
			radius * 2. + 4., Pars.camera_spider_width * radius );

		/* the NICMOS 2 spiders are slightly non-perpendicular to each other */

		if ( Pars.chip == NICMOS_2 || Pars.chip == NEWNICMOS_2 )
		{
			newimage = Rotate_image( temp1, dim, dim, 
					nicmos_x_offset-0.5, nicmos_y_offset-0.5, -0.6 );
			Free_image( temp1 );

			for ( y = 0; y < dim; ++y )
				for ( x = 0; x < dim; ++x )
					temp[y][x] *= newimage[y][x];

			Free_image( newimage );
		}
		else
		{
			for ( y = 0; y < dim; ++y )
				for ( x = 0; x < dim; ++x )
				temp[y][x] *= temp1[y][x];

			Free_image( temp1 );
		}

         	Dark_rectangle( temp, dim, dim, nicmos_x_offset, nicmos_y_offset, 
			Pars.camera_spider_width * radius, radius * 2. + 4. );

        	if ( Pars.chip == NICMOS_2 || Pars.chip == NICMOS_3 ||
		     Pars.chip == NEWNICMOS_2 || Pars.chip == NEWNICMOS_3 )
	        {
         		/* Paint in the pad covers  */

	         	for (i = 0; i < 3; i++ )
        	    		Dark_rot_rect( temp, dim, dim,
                	         	Pars.nic_pad_v3[i] * radius + nicmos_x_offset,
                        	 	Pars.nic_pad_v2[i] * radius + nicmos_y_offset,
                         		Pars.nic_pad_x[i] * radius,
                         		Pars.nic_pad_y[i] * radius,
	                         	Pars.nic_pad_rot[i] * M_PI / 180.0 );
        	}

		/* Rotate mask, if necessary, and multiply OTA pupil */

		newimage = Rotate_image( temp, dim, dim, 
					nicmos_x_offset-0.5, nicmos_y_offset-0.5, 
					Pars.nicmos_mask_rotation-Pars.ota_offset );
		Free_image( temp );
		for ( y = 0; y < dim; ++y )
			for ( x = 0; x < dim; ++x )
				pupil[y][x] *= newimage[y][x];
		Free_image( newimage );
	}

        /* WFC3 IR channel cold mask obscuration - note no pad covers */

	if ( Pars.chip == WFC3_IR )
        {
		wfc3_x_offset = Pars.wfc3_mask_x * radius + center;
		wfc3_y_offset = Pars.wfc3_mask_y * radius + center;

         	Dark_outside_circle( pupil, dim, dim, wfc3_x_offset, wfc3_y_offset, 
			Pars.camera_entrance_radius * radius );
         	Dark_inside_circle( pupil, dim, dim, wfc3_x_offset, wfc3_y_offset, 
			Pars.camera_secondary_radius * radius );
         	Dark_rectangle( pupil, dim, dim, wfc3_x_offset, wfc3_y_offset, 
			radius * 2. + 4., Pars.camera_spider_width * radius );
         	Dark_rectangle( pupil, dim, dim, wfc3_x_offset, wfc3_y_offset, 
			Pars.camera_spider_width * radius, radius * 2. + 4. );
	}

	/* Rotate the pupil */

	if ( Pars.theta != 0.0 )
	{
		newimage = Rotate_image( pupil, dim, dim,
				center-0.5, center-0.5, Pars.theta+Pars.ota_offset );

		for ( y = 0; y < dim; ++y )
			for ( x = 0; x < dim; ++x )
				pupil[y][x] = newimage[y][x];

		Free_image( newimage );
	}

	/* Add in extra WF/PC-1 WFC vignetting */

	if ( Pars.chip == 1 || Pars.chip == 2 || Pars.chip == 4 )
		Wfc_vignette( pupil, dim, dim );

	/* STIS CCD Lyot Stop */

	if ( Pars.chip == STIS_CCD )
	{
		newimage = Alloc_image( dim, dim );
		stis_x_offset = Pars.stis_stop_x_offset * radius + center;
		stis_y_offset = Pars.stis_stop_y_offset * radius + center;
		Draw_ellipse( newimage, dim, dim, stis_x_offset, stis_y_offset,
			Pars.stis_stop_x_radius*radius, Pars.stis_stop_y_radius*radius );
		Smooth( newimage, dim, dim, dim/512*5 );
		for ( y = 0; y < dim; ++y )
		   for ( x = 0; x < dim; ++x )
			pupil[y][x] *= newimage[y][x];
		Free_image( newimage );
	}

} /* Compute_pupil */


