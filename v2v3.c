/*  File     :  v2v3.c
 *  
 *  Contents :
 *      XY_to_V2V3 : Convert detector X,Y to pupil V2,V3; used by v2v3()
 *	V2V3_to_XY : Convert field V2,V3 to detector X,Y; used by v2v3()
 *	v2v3 : Compute pupil V2,V3 for a given detector X,Y
 *
 *  Author   :  John Krist 
 *  Date     :  October 2000
 *
 * Modified to include WFC3, Richard Hook & Felix Stoehr, March 2008
 *
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/*------------------------------------------------------------------------
*  XY_to_V2V3 :
*
*	Convert geometrically distorted detector X,Y pixel position to 
*	V2,V3 pupil coordinates.
*
*  Inputs :
*	x, y : Detector pixel position : (0,0) for lower left corner
*     camera : camera ID
*  Outputs :
*	v2, v3 : Position in arcsec in pupil coordinate system
*
*  NOTE : This routine only works for the ACS & WFC3
*-------------------------------------------------------------------------*/
void XY_to_V2V3( float x, float y, float *v2, float *v3, int camera )
{
	float	xc, yc, x2, x3, x4, y2, y3, y4, root2;

        root2 = sqrt(2.0);

	/* transforms are relative to the reference position */

	x = x - Pars.x_ref;
	y = y - Pars.y_ref;

	x2 = x * x;
	x3 = x2 * x;
	x4 = x3 * x;
	y2 = y * y;
	y3 = y2 * y;
	y4 = y3 * y;

	/* transform distorted x,y to corrected x,y (offset in arcsec from reference point) */

	xc  =  y*Pars.xy_to_xc[0] + x*Pars.xy_to_xc[1] +
	       y2*Pars.xy_to_xc[2] + x*y*Pars.xy_to_xc[3]  + x2*Pars.xy_to_xc[4] +
	       y3*Pars.xy_to_xc[5] + x*y2*Pars.xy_to_xc[6]  + x2*y*Pars.xy_to_xc[7] + x3*Pars.xy_to_xc[8] +
	       y4*Pars.xy_to_xc[9] + x*y3*Pars.xy_to_xc[10] + x2*y2*Pars.xy_to_xc[11] + x3*y*Pars.xy_to_xc[12] + x4*Pars.xy_to_xc[13];

	yc  =  y*Pars.xy_to_yc[0] + x*Pars.xy_to_yc[1] +
	       y2*Pars.xy_to_yc[2] + x*y*Pars.xy_to_yc[3]  + x2*Pars.xy_to_yc[4] +
	       y3*Pars.xy_to_yc[5] + x*y2*Pars.xy_to_yc[6]  + x2*y*Pars.xy_to_yc[7] + x3*Pars.xy_to_yc[8] +
	       y4*Pars.xy_to_yc[9] + x*y3*Pars.xy_to_yc[10] + x2*y2*Pars.xy_to_yc[11] + x3*y*Pars.xy_to_yc[12] + x4*Pars.xy_to_yc[13];

	/* in HRC/SBC, field +V2 is left, +V3 is up; for WFC, +V2 is right, +V3 is down. */
	/* Pars.v2 and Pars.v3 are in pupil coordinate system. */

        /* For WFC3 there is a 45 degree rotation as well as a flip */

	if ( camera == ACS_HRC || camera == ACS_HRC_OFFSPOT || camera == ACS_SBC )
	{
		*v2 = -Pars.v2 - xc;
		*v3 = -Pars.v3 + yc;
	}
        else if ( camera == WFC3_IR || camera == WFC3_UVIS1 || camera == WFC3_UVIS2 )
        {
                *v2 = -Pars.v2 + (-xc + yc)/root2;
                *v3 = -Pars.v3 + (xc + yc)/root2;
        }
	else
	{
		*v2 = -Pars.v2 + xc;
		*v3 = -Pars.v3 - yc;
	}
  
	/* convert to V2,V3 pupil coordinate system */

	*v2 = -(*v2);
        *v3 = -(*v3);

        /*
        printf("X,Y, V2,V3 (pupil): %f %f %f %f\n",x,y,*v2,*v3);
        */
}

/*------------------------------------------------------------------------
*  V2V3_to_XY :
*
*	Convert V2,V3 field coordinates to geometrically distorted detector
*	X,Y pixel position 
*
*  Inputs :
*	v2, v3 : Field position in arcsec (NOT in pupil coordinate system)
*       camera : Camera ID
*  Outputs :
*	x, y : Detector pixel position
*
*  NOTE : This routine only works for the ACS & WFC3
*-------------------------------------------------------------------------*/
void V2V3_to_XY( float *x, float *y, float v2, float v3, int camera )
{
	float	xc, yc, xc2, yc2, xc3, yc3, xc4, yc4, root2;

       root2 = sqrt(2.0);

	/* Pars.v2 and Pars.v3 are in pupil (not field coordinates) */

	if ( camera == ACS_HRC || camera == ACS_HRC_OFFSPOT || camera == ACS_SBC )
	{
		xc = (-Pars.v2) - v2;
		yc = v3 - (-Pars.v3);
	}
        else if ( camera == WFC3_IR || camera == WFC3_UVIS1 || camera == WFC3_UVIS2 )
        {
                xc = -Pars.v2 + (-v2 + v3)/root2;
                yc = -Pars.v3 + (v2 + v3)/root2;
        }
	else
	{
		xc = v2 - (-Pars.v2);
		yc = (-Pars.v3) - v3;
	}

	xc2 = xc * xc;
	xc3 = xc2 * xc;
	xc4 = xc3 * xc;
	yc2 = yc * yc;
	yc3 = yc2 * yc;
	yc4 = yc3 * yc;

	*x =    yc*Pars.xcyc_to_x[0] + xc*Pars.xcyc_to_x[1] +
	     	yc2*Pars.xcyc_to_x[2] + xc*yc*Pars.xcyc_to_x[3] + xc2*Pars.xcyc_to_x[4] + 
		yc3*Pars.xcyc_to_x[5] + xc*yc2*Pars.xcyc_to_x[6] + xc2*yc*Pars.xcyc_to_x[7] + xc3*Pars.xcyc_to_x[8] +
		yc4*Pars.xcyc_to_x[9] + xc*yc3*Pars.xcyc_to_x[10] + xc2*yc2*Pars.xcyc_to_x[11] + xc3*yc*Pars.xcyc_to_x[12] +
		xc4*Pars.xcyc_to_x[13];

	*y =    yc*Pars.xcyc_to_y[0] + xc*Pars.xcyc_to_y[1] +
	     	yc2*Pars.xcyc_to_y[2] + xc*yc*Pars.xcyc_to_y[3] + xc2*Pars.xcyc_to_y[4] + 
		yc3*Pars.xcyc_to_y[5] + xc*yc2*Pars.xcyc_to_y[6] + xc2*yc*Pars.xcyc_to_y[7] + xc3*Pars.xcyc_to_y[8] +
		yc4*Pars.xcyc_to_y[9] + xc*yc3*Pars.xcyc_to_y[10] + xc2*yc2*Pars.xcyc_to_y[11] + xc3*yc*Pars.xcyc_to_y[12] +
		xc4*Pars.xcyc_to_y[13];

	/* transforms are relative to reference position */

	*x = *x + Pars.x_ref;
	*y = *y + Pars.y_ref;
}

/*-------------------------------------------------------------------------
*  v2v3 :
*	Compute V2,V3 pupil coordinates for a given detector X,Y position;
*	pupil coordinates are the negative of field coordinates
*
*  Inputs :
*	camera : camera id number
*	  x, y : position on detector 
*  Outputs :
*	v2, v3 : v2,v3 pupil coordinates
*-------------------------------------------------------------------------*/
void v2v3( int camera, int x, int y, float *v2, float *v3 )
{
	float	t;

	if ( camera >= 1 && camera <= 8 )
	{
		/* WFPC1 coordinates are 0-799 */

		/* convert to arcsec from detector center */ 

		x = (x - 412) * Pars.pixel_size;	
		y = (y - 412) * Pars.pixel_size;

		/* adjust for detector rotation with respect to OTA */

		t = (90 + Pars.theta) * M_PI / 180.0;
		*v2 = -(x * cos(t) - y * sin(t)) + Pars.v2;
		*v3 = x * sin(t) + y * cos(t) + Pars.v3;
	}
	else if ( camera >= WFPC2_PC && camera <= WFPC2_WFC_4 )
	{
		/* WFPC2 coordinates are 0-799 */

		/* convert to arcsec from detector center */ 

		x = (x - 420) * Pars.pixel_size;	
		y = (y - 420) * Pars.pixel_size;

		/* adjust for detector rotation with respect to OTA */

		t = (90 + Pars.theta) * M_PI / 180.0;
		*v2 = -(x * cos(t) - y * sin(t)) + Pars.v2;
		*v3 = x * sin(t) + y * cos(t) + Pars.v3;
	}
	else if ( (camera >= NICMOS_1 && camera <= NICMOS_3) || (camera >= NEWNICMOS_1 && camera <= NEWNICMOS_3) )
	{
		/* NICMOS coordinates are 0-255 */

		/* convert to arcsec from detector center */ 

		x = (x - 128) * Pars.pixel_size;	
		y = (y - 128) * Pars.pixel_size;

		/* adjust for detector rotation with respect to OTA */

		t = (90 + Pars.theta) * M_PI / 180.0;
		*v2 = -(x * cos(t) - y * sin(t)) + Pars.v2;
		*v3 = x * sin(t) + y * cos(t) + Pars.v3;
	}
	else if ( camera >= ACS_WFC1 && camera <= WFC3_IR )
		XY_to_V2V3( (float)x, (float)y, v2, v3, camera );
	else
	{
		*v2 = Pars.v2;
		*v3 = Pars.v3;
	}

} /* v2v3 */

