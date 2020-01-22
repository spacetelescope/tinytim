/*  File     :  opd.c
 *  
 *  Contents :
 *            Compute_opd : Compute the OPD array for 547 nm.
 *
 *  Author   :  John Krist (STScI)
 *  Date     :  January 1992
 *
 *  Modifications :
 *      March 1994 - JEK
 *        Modified ComputeOPD to adjust focus depending on PSF position for
 *        WFPC2 cameras only.
 *
 *      November 1996 - JEK
 * 	  Removed support for old mirror maps.
 *
 *      June 1997 - JEK
 *        Added ComputeNicmosAber() which computes the field- and focus-
 *        dependent NICMOS aberrations.
 *
 *	December 1998 - JEK
 *	  Added support for WFPC2 field dependent aberrations (Field_aberration)
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*--------------------------------------------------------------------------
*  Compute_opd  :
*       Compute a 2-d image which contains the optical path differences
*    (given in waves at 547 nm) for each point in the pupil. For each point,
*   the OPD is computed as the sum of the Zernike polynomials.  For a table
*    of the Zernike polynomials, see the HST OTA Handbook (STScI).
*       Since the pupil is defined only within the central N/2 by N/2 portion
*    of the image, only those points within that region are computed.
*       If a WFPC2 camera is modelled, then the focus is adjusted based on
*    PSF position to account for change in focus across each chip.
*
*    Inputs:
*       dim   :   Dimension of one side of Opd array.
*       opd   :   Image array allocated with Alloc_image of size Dim by Dim.
*       psf_x,
*       psf_y  :   Position of PSF on detector (only used for WFPC2).
*
*   Returns:
*       Returns the OPD function in waves at 547 nm in Opd.
*
*   Changes :  Adjusted distance computation by 0.5 pixels for correct results,
*               but this does not change things very much.
*
*  The x,y slopes for WFPC2 are from the WFPC2 Science Calibration Report.
*  The conversion from encoder steps/pixel to waves/pixel @ 547 nm used was :
*       Dmm = SlopeEnc * 0.000467 * 25.4 = Slope in f24 focal plane in mm/pix.
*       SlopeMicrons = (1000. * Dmm) / (3.887*8*24^2) = Zernike microns/pix
*       SlopeWave = SlopeMicrons / 0.547 = Zernike waves/pix @ 547 nm.
*-----------------------------------------------------------------------------*/
void Compute_opd( int dim, float **opd, int psf_x, int psf_y )
{
	float   center, radius, v3dist[4096], v3sqr[4096], v2sqr;
	float   r, r2, r3, r4, r5, r6, theta, value, v2dist;
	float   sin_theta, sin_theta2, sin_theta3, sin_theta4, sin_theta5;
	float   cos_theta, cos_theta2, cos_theta3, cos_theta4, cos_theta5;
	float   *zval, **map, **temp_opd, v2_offset, v3_offset;
	int     chip, x, y, x1, x2, y1, y2, v2, v3;


	Compute_aberrations( psf_x, psf_y );

	zval = Pars.zval;

	center = dim / 2.0;
	radius = dim / 4.0;

	/* Since the pupil exists only in the central Dim/2 region,  *
	 * compute only within this region.                          */

	x1 = y1 = (int)(center - radius) - 6;
	x2 = y2 = (int)(center + radius) + 6;

	/* Compute lookup table of distances */

	for ( x = 0; x < dim; ++x )
	{
		v3dist[x] = (x - center) / radius;
		v3sqr[x] = v3dist[x] * v3dist[x];
	}

	for ( v2 = y1; v2 <= y2; ++v2 )
	{
		v2dist = v3dist[v2];
		v2sqr = v3sqr[v2];

		for ( v3 = x1; v3 <= x2; ++v3 )
		{
			value = 0.0;

			r = sqrt( v3sqr[v3] + v2sqr );
			r2 = r * r;
			r3 = r2 * r;
			r4 = r3 * r;
			r5 = r4 * r;
			r6 = r5 * r;

			if ( v3dist[v3] != 0.0 || v2dist != 0.0 )
				theta = atan2( v3dist[v3], v2dist );
			else
				theta = 0.0;

			sin_theta = sin(theta);
			cos_theta = cos(theta);
			sin_theta2 = sin(theta*2.);
			cos_theta2 = cos(theta*2.);
			sin_theta3 = sin(theta*3.);
			cos_theta3 = cos(theta*3.); 
			sin_theta4 = sin(theta*4.);
			cos_theta4 = cos(theta*4.);
			sin_theta5 = sin(theta*5.);
			cos_theta5 = cos(theta*5.);

			if ( zval[2] != 0.0 )    /* X tilt */
				value += zval[2] * 1.8992573 * r * cos_theta;
			if ( zval[3] != 0.0 )    /* Y tilt */
				value += zval[3] * 1.8992573 * r * sin_theta;
			if ( zval[4] != 0.0 )    /* focus */
				value += Pars.focus * 3.8874443 * 
					   (r2 - 0.554450);
			if ( zval[5] != 0.0 )    /* 0 degree astigmatism */
				value += Pars.xastig * 2.3137662 * r2 * 
					   cos_theta2;
			if ( zval[6] != 0.0 )    /* 45 degree astigmatism */
				value += Pars.yastig * 2.3137662 * r2 * 
					   sin_theta2;
			if ( zval[7] != 0.0 )    /* X coma */
				value += Pars.xcoma * 8.3345629 * 
					(r3 - 0.673796 * r) * cos_theta;
			if ( zval[8] != 0.0 )    /* Y coma */
				value += Pars.ycoma * 8.3345629 * 
					(r3 - 0.673796 * r) * sin_theta;
			if ( zval[9] != 0.0 )    /* X clover */
				value += Pars.xclover * 2.6701691 * r3 * cos_theta3;
			if ( zval[10] != 0.0 )   /* Y clover */
				value += Pars.yclover * 2.6701691 * r3 * sin_theta3;
			if ( zval[11] != 0.0 )   /* 3rd order spherical */
				value += Pars.spherical * 16.895979 * 
					(r4 - 1.1089 * r2 + 0.241243);
			if ( zval[12] != 0.0 )   /* spherical astigmatism */
				value += zval[12] * 12.033645 * 
					(r4 - 0.750864 * r2) * cos_theta2;
			if ( zval[13] != 0.0 )   /* 45 degree spherical astigmatism */
				value += zval[13] * 12.033645 * 
					(r4 - 0.750864 * r2) * sin_theta2;
			if ( zval[14] != 0.0 )   /* Ashtray */
				value += zval[14] * 2.9851527 * r4 * cos_theta4;
			if ( zval[15] != 0.0 )   /* Ashtray */
				value += zval[15] * 2.9851527 * r4 * sin_theta4;
			if ( zval[16] != 0.0 )
				value += zval[16] * 36.321412 * (r5 - 1.230566 
					* r3 + 0.323221 * r) * cos_theta;
			if ( zval[17] != 0.0 )
				value += zval[17] * 36.321412 * (r5 - 1.230566 
					* r3 + 0.323221 * r) * sin_theta;
			if ( zval[18] != 0.0 )
				value += zval[18] * 16.372202 * (r5 - 0.8001 
					* r3) * cos_theta3;
			if ( zval[19] != 0.0 )
				value += zval[19] * 16.372202 * (r5 - 0.8001 
					* r3) * sin_theta3;
			if ( zval[20] != 0.0 )
				value += zval[20] * 3.2700486 * r5 * cos_theta5;
			if ( zval[21] != 0.0 )
				value += zval[21] * 3.2700486 * r5 * sin_theta5;
			if ( zval[22] != 0.0 )   /* 5th order spherical */
				value += zval[22] * 74.82446 * (r6 - 1.663350 
					* r4 + 0.803136 * r2 - 0.104406);

			opd[v2][v3] = value;
		}
	}

	chip = Pars.chip;

	/* Rotate the OPD to match the instrumental rotation */

	if ( Pars.theta != 0.0 )
	{
		temp_opd = Rotate( opd, dim, dim, center, center, 
				dim, dim, center, center, Pars.theta, 1.0 );
		for ( y = 0; y < dim; ++y )
			for ( x = 0; x < dim; ++x )
				opd[y][x] = temp_opd[y][x];
		Free_image( temp_opd );
	}

	/* Add in mirror maps.  The mirror map routines return rotated maps. */

	if ( Pars.use_map != 0 )	/* Use maps */
	{
		v2v3( chip, psf_x, psf_y, &v2_offset, &v3_offset );
		map = Compute_map( dim, chip, Pars.theta, v2_offset, v3_offset );

		for ( y = y1; y <= y2; ++y )
			for ( x = x1; x <=x2; ++x )
				opd[y][x] += map[y][x];

		Free_image( map );
	}

} /* Compute_opd */

