/* File      : zernike.c
 *
 * Contents  :
 * 	Julian_day : Convert date to Julian day
 *    Compute_nicmos_aber : Compute Nicmos focus & field dependent aberrations
 *       Field_aberration : Compute field dependent aberration
 *   Wfc3_field_aberration : Compute WFC3 field dependent aberration
 *   Acs_field_aberration : Compute ACS field dependent aberration
 *    Compute_aberrations : Compute field/focus dependent aberrations
 *	Nicmos_aberrations : Determine aberrations for NICMOS
 *      Set_old_aberrations : Set the aberrations for WF/PC-1 and pre-COSTAR FOC.
 *      Read_zernike_data  : Read Zernike coefficients from a file.
 *      Get_aberrations : Determine aberrations.
 *
 * Author    : John Krist
 * Date      : May 1993
 *
 * Modifications :
 *
 *      August 1993 - JEK
 *        Added FOC focus change on Oct 7, 1992, to SetAberrations.
 *
 *      March 1994 - JEK
 *        Modified ReadZernikeTable to account for separate WFPC2 tables.
 *
 *      February 1996 - Richard Hook
 *        Added NICMOS
 *
 *	Dec 1997 - JEK
 *	  Added STIS, ACS
 *
 *      May 1999 - JEK
 *	  Added NicmosAberrations
 *
 *      March 2008 - RNH/FS
 *        Added WFC3 support
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "tinytim.h"

#define NUM_DATES  20

static char *Aberration_names[] = { "Z1 : Constant", "Z2 : X Tilt", 
	"Z3 : Y Tilt", "z4 : Focus", "z5 : 0 deg Astigmatism", 
	"z6 : 45 deg Astigmatism", "z7 : V2 Coma", "z8 : V3 Coma",
	"z9 : X Clover", "z10 : Y Clover", "z11 : 3rd order Spherical",
	"Z12 : X Spher. Astig", "Z13 : Y Spher. Astig",
	"Z14 : X Ashtray", "Z15 : Y Ashtray", "Z16", "Z17", "Z18",
	"Z19", "Z20", "Z21", "z22 : 5th order Spherical" };

/*------------------------------------------------------------------------------
*  Julian_day :
*	Return Julian day of specified date
*
*  Inputs :
*	day, month, year : Date in numeric format
*
*  Returns :
*	Julian day.
*---------------------------------------------------------------------------*/
static long Julian_day( int day, int month, int year, int print_date )
{
        long    greg, jy, jm, ja, jul;
	static  char *months[12] = {"January", "February", "March", "April",
			"May", "June", "July", "August", "September", 
			"October", "November", "December" };

        greg = 15 + 31 * (10 + 12 * (long)1582);
        if ( month > 2 )
        {
                jy = year;
                jm = month + 1;
        }
        else
        {
                jy = year - 1;
                jm = month + 13;
        }

        jul = (long)(365.25 * jy) + (long)(30.6001 * jm) + day + 1720995;
        if ( (day + 31 * (month + 12 * year)) >= greg )
        {
                ja = (long)(0.01 * jy);
                jul = jul + 2 - ja + (long)(0.25 * ja);
        }

	if ( print_date == 1 )
		printf("\nDate is %s %d, %d\n", months[month-1], day, year ); 

        return( jul );
}

/*-------------------------------------------------------------------------
*  Compute_nicmos_aber :
*
*  Determine the Zernike coefficients for the field and focus dependent
*  aberrations.
*
*  Inputs :
*	psf_x, psf_y : Position on detector (0-255)
*  Outputs :
*       Returns in Pars.focus, etc, the Zernike coefficients at 547 nm 
*       at position (psf_x,psf_y), corrected for defocus, in the V2,V3 system.
*
*  NOTE : The input values Pars.xastig, etc, are the aberrations in waves 
*         @ 547 nm at best focus (z4=0) for the center of the field in 
*	  (V2,V3) coordinates.  Pars.focus should be 0.0 or contain a 
*	  focus offset due to breathing.
*------------------------------------------------------------------------*/
static void Compute_nicmos_aber( int psf_x, int psf_y )
{
	float z4, z5, z6, z7, z8, z11;
	float f, ff, xy, t;
	int   i, j;

	/* Adjust aberrations with respect to focus.  Aberrations are in   *
	 * detector coordinate system in microns. Zero offsets are applied *
	 * later.                                                          */

	z4 = Pars.focus * 0.547;	/* convert focus to microns RMS */
        z5 = z6 = z7 = z8 = 0.0;

	if ( Pars.adjust_for_focus != 0 )
	{
		f = z4;
		ff = f * f;

		z5 = Pars.z5_vs_z4[0] * ff + Pars.z5_vs_z4[1] * f;  
		z6 = Pars.z6_vs_z4[0] * ff + Pars.z6_vs_z4[1] * f;  
		z7 = Pars.z7_vs_z4[0] * ff + Pars.z7_vs_z4[1] * f;  
		z8 = Pars.z8_vs_z4[0] * ff + Pars.z8_vs_z4[1] * f;  
		z11 = Pars.z11_vs_z4[0] * ff + Pars.z11_vs_z4[1] * f; 
	}
	else
	{
		z11 = 0.0;
	}

	/* Adjust aberrations with respect to field position.  Aberrations *
      	 * are in detector coordinate system in microns.                   */

	if ( Pars.adjust_for_xy != 0 )
	{
		for ( i = 0; i <= 2; ++i )
		{
		   for ( j = 0; j <= 2; ++j )
		   {
			xy = pow((float)psf_x, (float)i) * pow((float)psf_y, (float)j);
			z4 += Pars.z4_vs_xy[i][j] * xy; 
			z5 += Pars.z5_vs_xy[i][j] * xy; 
			z6 += Pars.z6_vs_xy[i][j] * xy; 
			z7 += Pars.z7_vs_xy[i][j] * xy; 
			z8 += Pars.z8_vs_xy[i][j] * xy; 
		   }
		}
	}

	/* Rotate aberrations to V2,V3 coordinate system and add    *
	 * center-of-field offsets (from the input file), converted *
	 * to waves @ 547 nm.					    */

	if ( (Pars.adjust_for_focus != 0) || (Pars.adjust_for_xy != 0) )
	{
		t = -Pars.theta * M_PI / 180.0;
		/* t = -M_PI / 4.0; */
		Pars.xastig += (z5 * cos(2.0*t) - z6 * sin(2.0*t)) / 0.547;
		Pars.yastig += (-z5 * sin(2.0*t) - z6 * cos(2.0*t)) / 0.547;
		Pars.xcoma  += (z7 * cos(t) - z8 * sin(t)) / 0.547;
		Pars.ycoma  += (-z7 * sin(t) - z8 * cos(t)) / 0.547;
	}

	Pars.focus = z4 / 0.547;       /* waves @ 547 nm */
	Pars.spherical = Pars.spherical + z11 / 0.547;

} /* Compute_nicmos_aber */

/*--------------------------------------------------------------------------
*  Field_aberration
*	Compute the specified, field dependent aberration for a given
*  	position on the detector, using quadratic relations.
*
*---------------------------------------------------------------------------*/
static float Field_aberration( int x, int y, int aberration )
{
	float *caber = &Pars.caber[aberration][0];

	return( caber[0]*x + caber[1]*y + caber[2]*x*x + caber[3]*y*y + caber[4]*x*y );
}

/*--------------------------------------------------------------------------
*  Wfc3_field_aberration
*	Compute the field dependent aberration for a given
*       position on the detector.
*
*---------------------------------------------------------------------------*/
static float Wfc3_field_aberration( int xdet, int ydet, int znum )
{

	int 	i, j;
	float	d, val, x, y, convert;

	val = 0.0;

	if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 ) {
		d       = 40.96;
		convert = 0.547/0.6328;
	} else {
	
		d       = 10.24;
		convert = 0.547;
	}	

	/* focus (znum=0) is defined relative to chip, not global, (x,y) system */

	if ( Pars.chip == WFC3_UVIS1) {
	       /* This was for the focus of the ACS camera. We do not have that, i.e.
	          we have given the focii in the same way as all the other coefficients
	          && znum != 0 )
	       */
	   ydet = ydet + 2048;
	}	

	x = xdet / d;      /* + 4; it is not clear where that came from in the ACS code */
	y = ydet / d;      /* + 4; it is not clear where that came from in the ACS code */

	/* compute field-dependent aberrations offsets in V2,V3 pupil system; *
         * coefficients are applied to detector x,y pixel coordinates and     *
	 * return aberrations in V2,V3 pupil system.  The field-              *
	 * dependent coefficients are in microns and have to be converted to  *
	 * waves at 547 nm.						      */
	
	for ( i = 0; i <= 5; ++i )
		for ( j = 0; j <= 5; ++j )
			val += Pars.wfc3_aber[znum][i][j] * pow((float)x,(float)i) * pow((float)y,(float)j);

        /* the aberrations are given in 628.3 nm for the UVIS chip and in 1 micron for the IR chip */	
	return( val / convert );

} /* Wfc3_field_aberration */

/*--------------------------------------------------------------------------
*  Acs_field_aberration
*	Compute the field dependent aberration for a given
*       position on the detector.
*
*---------------------------------------------------------------------------*/
static float Acs_field_aberration( int xdet, int ydet, int znum )
{
	int 	i, j;
	float	d, val, x, y;

	val = 0.0;

	if ( Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 )
		d = 40.96;
	else
		d = 10.24;

	/* focus (znum=0) is defined relative to chip, not global, (x,y) system */

	if ( Pars.chip == ACS_WFC1 && znum != 0 )
		ydet = ydet + 2048;

	x = xdet / d + 4;
	y = ydet / d + 4;

	/* compute field-dependent aberrations offsets in V2,V3 pupil system; *
         * coefficients are applied to detector x,y pixel coordinates and     *
	 * return aberrations in V2,V3 pupil system.  The field-              *
	 * dependent coefficients are in microns and have to be converted to  *
	 * waves at 547 nm.						      */

	for ( i = 0; i <= 5; ++i )
		for ( j = 0; j <= 5; ++j )
			val += Pars.acs_aber[znum][i][j] * pow((float)x,(float)i) * pow((float)y,(float)j);
	return( val / 0.547 );

} /* Acs_field_aberration */

/*--------------------------------------------------------------------------
* Compute_aberrations :
*	Compute field/focus dependent aberrations (for NICMOS, WFPC2,
*	or ACS.
*
* Inputs :
*	psf_x, psf_y : PSF position on detector in pixels
*-------------------------------------------------------------------------*/
void Compute_aberrations( int psf_x, int psf_y )
{
	float	xx, yy;

	/* Set focus/position dependent aberrations */

	/* get aberrations for field center */

	Pars.focus = Pars.zval[4];
	Pars.xastig = Pars.zval[5];
	Pars.yastig = Pars.zval[6];
	Pars.xcoma = Pars.zval[7];
	Pars.ycoma = Pars.zval[8];
	Pars.xclover = Pars.zval[9];
	Pars.yclover = Pars.zval[10];
	Pars.spherical = Pars.zval[11];

	/* add in offsets for field dependence */

	if ( Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4 && Pars.adjust_for_xy != 0 )
	{
		xx = psf_x - 400;
		yy = psf_y - 400;
		Pars.focus  += Field_aberration( xx, yy, 0 );
		Pars.xastig += Field_aberration( xx, yy, 1 );
		Pars.yastig += Field_aberration( xx, yy, 2 );
		Pars.xcoma  += Field_aberration( xx, yy, 3 );
		Pars.ycoma  += Field_aberration( xx, yy, 4 );
	}
	else if ( Pars.chip >= ACS_WFC1 && Pars.chip <= ACS_SBC && Pars.adjust_for_xy != 0 )
	{
		Pars.focus  += Acs_field_aberration( psf_x, psf_y, 0 );
		Pars.xastig += Acs_field_aberration( psf_x, psf_y, 1 );
		Pars.yastig += Acs_field_aberration( psf_x, psf_y, 2 );
		Pars.xcoma  += Acs_field_aberration( psf_x, psf_y, 3 );
		Pars.ycoma  += Acs_field_aberration( psf_x, psf_y, 4 );
	}
	else if ( Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3 )
	{
		Compute_nicmos_aber( psf_x, psf_y );
	}
	else if (Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 || Pars.chip == WFC3_IR)
	{
            	Pars.focus  += Wfc3_field_aberration( psf_x, psf_y, 0 );
		Pars.xastig += Wfc3_field_aberration( psf_x, psf_y, 1 );
		Pars.yastig += Wfc3_field_aberration( psf_x, psf_y, 2 );
		Pars.xcoma  += Wfc3_field_aberration( psf_x, psf_y, 3 );
		Pars.ycoma  += Wfc3_field_aberration( psf_x, psf_y, 4 );
	}
	else
	{
	   printf("No field dependent aberrations.");
	}

} /* Compute_aberrations */

/*------------------------------------------------------------------------------
*  Nicmos_aberrations :
*
*	Determine date-dependent NICMOS aberrations for field center, and
*	determine the cold mask offset.  This routine asks the user for 
*	the observation date.
*
*-----------------------------------------------------------------------------*/
static void Nicmos_aberrations( FILE *file, struct DateStruct *obs_date )
{
	float	*zval;
	float	z4_1, z5_1, z6_1, z7_1, z8_1, z9_1, z10_1, z11_1;
	float	z4_2, z5_2, z6_2, z7_2, z8_2, z9_2, z10_2, z11_2;
	float	xmask_1, ymask_1, xmask_2, ymask_2, pix_size_1, pix_size_2;
	float	mask_rot_1, mask_rot_2, ota_rot_1, ota_rot_2;
	float	diff, day_obs, day1, day2;
	int	date1, date2, month1, month2, year1, year2, used_pam;
	int	i, status;
	char	line[MAX_STRING];


	zval = &Pars.zval[0];

	printf("\nEnter the date of observation in the form dd mm yyyy.\n");
	printf("For example, 6 April 1998 would be \"6 4 1998\".\n");
	printf("Valid range is 15 April 1997 to 31 December 1998.\n");
	printf("To enter PAM positions instead, type PAM.\n\n");
	printf("Date : ");
	scanf( "%s", line );

	Pars.best_focus_pam = 0.0;
	Pars.observation_pam = 0.0;

	if ( strstr(line,"PAM") == NULL && strstr(line,"pam") == NULL )
	{
		used_pam = 0;
		sscanf(line, "%d", &obs_date->day);
		scanf("%d %d", &obs_date->month, &obs_date->year );
	}
	else
	{
		used_pam = 1;
		printf("Enter PAM position of best focus (mm) : ");
		scanf("%f", &Pars.best_focus_pam );
		printf("Enter PAM position of observation (mm) : ");
		scanf("%f", &Pars.observation_pam );

		printf("\nEnter the date of observation in the form dd mm yyyy.\n");
		printf("For example, 6 April 1998 would be \"6 4 1998\".\n");
		printf("Valid range is 15 April 1997 to 31 December 1998.\n");
		printf("This date will be used to determine the other aberrations.\n\n");
		printf("Date : ");
		scanf( "%d %d %d", &obs_date->day, &obs_date->month, &obs_date->year );
	}

	day_obs = Julian_day(obs_date->day,obs_date->month,obs_date->year,1);

	/* skip header lines in aberration table (begin with #) */

	line[0] = '#';
	while ( line[0] == '#' )
		fgets( line, MAX_STRING-2, file );
	
	/* read in time dependent aberrations (in microns) */

	sscanf( line, "%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",
	   	&date1, &month1, &year1, 
	   	&z4_1, &z5_1, &z6_1, &z7_1, &z8_1, &z9_1, &z10_1, &z11_1,
	   	&xmask_1, &ymask_1, &mask_rot_1, &ota_rot_1, &pix_size_1 );

	day1 = Julian_day(date1,month1,year1,0);

	if ( day_obs < day1 )
	{
		printf("ERROR : Valid range is 15 April 1997 to 31 December 1998.\n");
		fclose( file );
		exit(0);
	}

	status = fscanf( file, "%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",
			&date2, &month2, &year2, 
			&z4_2, &z5_2, &z6_2, &z7_2, &z8_2, &z9_2, &z10_2, &z11_2,
			&xmask_2, &ymask_2, &mask_rot_2, &ota_rot_2, &pix_size_2 );
	day2 = Julian_day(date2,month2,year2,0);

	while ( day_obs > day2 )
	{
		day1 = day2;
		z4_1 = z4_2;
		z5_1 = z5_2;
		z6_1 = z6_2;
		z7_1 = z7_2;
		z8_1 = z8_2;
		z9_1 = z9_2;
		z10_1 = z10_2;
		z11_1 = z11_2;
		xmask_1 = xmask_2;
		ymask_1 = ymask_2;
		mask_rot_1 = mask_rot_2;
		ota_rot_1 = ota_rot_2;
		pix_size_1 = pix_size_2;

		status = fscanf( file, "%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",
			   &date2, &month2, &year2, 
			   &z4_2, &z5_2, &z6_2, &z7_2, &z8_2, &z9_2, &z10_2, &z11_2,
			   &xmask_2, &ymask_2, &mask_rot_2, &ota_rot_2, &pix_size_2 );
		if ( status == EOF )
			break;
		day2 = Julian_day(date2,month2,year2,0);
	}

	if ( day_obs > day2 )
	{
		printf("ERROR : Valid range is 15 April 1997 to 31 December 1998.\n");
		fclose( file );
		exit(0);
	}

	/* Interpolate between dates to determine aberrations (for field center)  *
	 * and cold mask offsets.  Aberrations & mask shifts are in pupil V2,V3.  *
	 * Field dependent aberrations are computed during PSF generation, since  *
	 * PSFs at multiple positions can be generated in a single run.           */

	diff = day_obs - day1;

	/* all aberrations in microns */

	if ( used_pam )
		zval[4] = Pars.pam_to_z4 * (Pars.observation_pam - Pars.best_focus_pam);
	else
		zval[4] =  diff * (z4_2 - z4_1) / (day2 - day1) + z4_1;

	zval[5] =  diff * (z5_2 - z5_1) / (day2 - day1) + z5_1;
	zval[6] =  diff * (z6_2 - z6_1) / (day2 - day1) + z6_1;
	zval[7] =  diff * (z7_2 - z7_1) / (day2 - day1) + z7_1;
	zval[8] =  diff * (z8_2 - z8_1) / (day2 - day1) + z8_1;
	zval[9] =  diff * (z9_2 - z9_1) / (day2 - day1) + z9_1;
	zval[10] = diff * (z10_2 - z10_1) / (day2 - day1) + z10_1;
	zval[11] = diff * (z11_2 - z11_1) / (day2 - day1) + z11_1;
	Pars.nicmos_mask_x_off = diff * (xmask_2 - xmask_1) / (day2 - day1) + xmask_1;
	Pars.nicmos_mask_y_off = diff * (ymask_2 - ymask_1) / (day2 - day1) + ymask_1;
	Pars.nicmos_mask_rotation = diff * (mask_rot_2 - mask_rot_1) / (day2 - day1) + mask_rot_1;
	Pars.ota_offset = diff * (ota_rot_2 - ota_rot_1) / (day2 - day1) + ota_rot_1;
	Pars.pixel_size = diff * (pix_size_2 - pix_size_1) / (day2 - day1) + pix_size_1;

	/* Create coefficients table needed by ReadZernikeData; convert to waves @ 547 nm */

	strcpy( line, "#\n" );
	Store_entry( line );
	strcpy( line, "#  Zernike Polynomial Coefficients for field center\n" );
	Store_entry( line );
	strcpy( line, "#\n" );
	Store_entry( line );
	strcpy( line, " 547.0 # Reference wavelength (nm) for coeffs\n" );
	Store_entry( line );
	strcpy( line, "  22    # Last Zernike in table\n" );
	Store_entry( line );

	for ( i = 1; i <= LASTZERNIKE; ++i )
	{
		sprintf( line, "%f # %s\n", zval[i]/0.547, Aberration_names[i-1] );
		Store_entry( line );
	}

} /* Nicmos_aberrations */

/*------------------------------------------------------------------------------
*  Set_old_aberrations :
*
*	Determine aberrations for WF/PC-1 or pre-COSTAR FOC for a given date.
*	The aberrations for certain dates are read from a file and the values
*	for the observation date are computed, including mirror moves and
*	desorption.
*
*  Inputs/Outputs :
*       zval = Floating point array containing Zernike polynomial coefficients.
*               The Zernike values for the given date are placed in zval.
*       chip = chip/Camera code.
*	obs_date = Date of observation (returned).
*	table_name = Name of table containing aberration information.
*---------------------------------------------------------------------------*/ 
static void Set_old_aberrations( struct DateStruct *obs_date, FILE *file )
{
	int	focus_day, focus_month, focus_year;
	int     i, j, day[NUM_DATES];
	int     change_date, day_number;
	long	launch_date;
	float   z4[NUM_DATES], z5[NUM_DATES], z6[NUM_DATES], z7[NUM_DATES];
	float   z8[NUM_DATES], z9[NUM_DATES], z10[NUM_DATES], *zval;
	float	z11, z22, desorption, focus_offset, zero_z4, zero_sm;
	char    entry[MAX_STRING];

	zval = &Pars.zval[0];

	/* HST was launched on April 24, 1990 */

	launch_date = Julian_day( 24, 4, 1990, 0 );

	/* Get the date of the observation.  This will be used to select *
	 * the aberrations appropriate for the secondary mirror position *
	 * on that date and to compute defocus due to desorption.	 */

	printf("\nEnter the date of observation in the form dd mm yyyy.\n");
	printf("For example, 6 April 1991 would be \"6 4 1991\".\n");
	printf("Valid range is 24 April 1990 to 2 December 1993.\n\n");
	printf("Date : ");
	scanf("%d %d %d", &obs_date->day, &obs_date->month, &obs_date->year );

	day_number = Julian_day(obs_date->day,obs_date->month,obs_date->year,1) - launch_date;

	sscanf( Get_string(file), "%d", &j );   /* Number of entries */

	/* Get z4 and secondary mirror position on August 15, 1990 */

	sscanf( Get_string(file), "%f %f", &zero_z4, &zero_sm );

	/* Read in aberrations for each secondary mirror position */

	for ( i = 0; i < j; ++i )
	{
		sscanf( Get_string(file), 
			"%d %d %d %f %f %f %f %f %f %f %f %f",
			&focus_day, &focus_month, &focus_year,
			&z4[i], &z5[i], &z6[i], &z7[i], &z8[i], &z9[i], &z10[i],
			&z11, &z22 );
		day[i] = Julian_day(focus_day,focus_month,focus_year,0) - launch_date;
	}

	/* z4[i] is the absolute OTA secondary mirror position in microns. *
	 * The Zernike focus constant for August 15, 1990, is read in, as  *
	 * well as the mirror position on that date.  This date is the     *
	 * zero point for the desorption curve.	 The focus for a given     *
	 * date is determined by adding to the zero point focus the        *
	 * desorption which has occured from Aug. 15, 1990 to that date,   *
	 * then subtracting the intervening mirror adjustments.		   *
	 * Note that zero_z4 is the focus for the PC in all cases.	   *
	 * The other aberrations are Zernike coefficients for that date.   *
	 * z11 & z22 (3rd & 5th order spherical) are constant with time.   */

	for ( i = 1; i < j; ++i )
		if ( day_number < day[i] )
			break;

	i = i - 1;

	zval[5] = z5[i];	/* 0 degree astigmatism */
	zval[6] = z6[i];	/* 45 degree astigmatism */
	zval[7] = z7[i];	/* X (V2) Coma */
	zval[8] = z8[i];	/* Y (V3) Coma */
	zval[9] = z9[i];	/* X clover */
	zval[10] = z10[i];	/* Y clover */
	zval[11] = z11;		/* 3rd order Spherical */
	zval[22] = z22;		/* 5th order Spherical */

	/* desorption is the change in the distance of the secondary       *
	 * mirror to the primary in microns due to desorption, relative to *
	 * the separation on August 16, 1990.				   */

	desorption = 111.245 - 93.229 * exp(-day_number / 911.095) - 
		77.387 * exp(-day_number / 115.773);

	/* Adjust for secondary mirror moves */

	focus_offset = desorption - (z4[i] - zero_sm);

	/* The FOC was separately refocused on October 7, 1992 by the  *
	 * equivalent of 6.7 microns of OTA secondary mirror movement. *
	 * There was also an apparent 8.9 micron focus offset between  *
	 * the FOC and WFPC1.					       */

	if ( Pars.chip == FOC_F48 || Pars.chip == FOC_F96 )
	{
		focus_offset += 8.9;
		change_date = Julian_day(7, 10, 1992, 0 ) - launch_date;
		if ( day_number >= change_date )
			focus_offset += 6.7;
	}

	/* Convert microns of secondary mirror movement to change in the *
	 * Zernike focus coefficient in waves @ 547 nm.			 */

	focus_offset = focus_offset * 110.0 / (3.887 * 8.0 * 24.0 * 24.0 * 0.547);

	/* Add focus offset to zeropoint */

	zval[4] = zero_z4 - focus_offset;

	/* Create coefficients table needed by ReadZernikeData */

	strcpy( entry, "#\n" );
	Store_entry( entry );
	strcpy( entry, "#  Zernike Polynomial Coefficients\n" );
	Store_entry( entry );
	strcpy( entry, "#\n" );
	Store_entry( entry );
	strcpy( entry, " 547.0 # Reference wavelength (nm) for coeffs\n" );
	Store_entry( entry );
	strcpy( entry, "  22    # Last Zernike in table\n" );
	Store_entry( entry );

	for ( i = 1; i <= LASTZERNIKE; ++i )
	{
		sprintf( entry, "%f # %s\n", zval[i], Aberration_names[i-1] );
		Store_entry( entry );
	}

} /* Set_old_aberrations */

/*--------------------------------------------------------------------------
*  Read_zernike_data :
*       Read Zernike values from a file.  
*
*  Inputs :
*       file  :  Pointer to aberration file (already opened).
*--------------------------------------------------------------------------*/
void Read_zernike_data( FILE *file )
{
	int     i, last_zernike;
	float   ref_wave;

	sscanf( Get_entry(file), "%f", &ref_wave );

	ref_wave = ref_wave / 1000.0;     /* Convert from nm to microns */

	sscanf( Get_entry(file), "%d", &last_zernike );
	for ( i = 1; i <= last_zernike; ++i )
	{
		sscanf( Get_entry(file), "%f", &Pars.zval[i] );
		Pars.zval[i] *= ref_wave / 0.547;
	}

} /* Read_zernike_data */

/*--------------------------------------------------------------------------
*  Get_aberrations :
*
*       Open aberration file for a given camera and compute aberration
*	values.  Some cameras will require the user to enter the
*	observation date for time-dependent aberrations.  If so, the
*	routine will return the entered date in obs_date. 
*
*--------------------------------------------------------------------------*/
void Get_aberrations( struct DateStruct *obs_date )
{
	FILE    *file;
	int     i;
	char	table_name[MAX_STRING];

	for ( i = 0; i <= LASTZERNIKE; ++i )
		Pars.zval[i] = 0.0;

	/* Default_dir adds the Tiny Tim directory path to the filename */

	switch ( Pars.chip )
	{
		     case 1 :   Default_dir("wfpcwf1.tab", table_name);  break;
		     case 2 :   Default_dir("wfpcwf2.tab", table_name);  break;
		     case 3 :   Default_dir("wfpcwf3.tab", table_name);  break;
		     case 4 :   Default_dir("wfpcwf4.tab", table_name);  break;
	             case 5 :   Default_dir("wfpcpc5.tab", table_name);  break;
		     case 6 :   Default_dir("wfpcpc6.tab", table_name);  break;
		     case 7 :   Default_dir("wfpcpc7.tab", table_name);  break;
		     case 8 :   Default_dir("wfpcpc8.tab", table_name);  break;
	      case WFPC2_PC :   Default_dir("wfpc2pc.tab", table_name);  break;
	   case WFPC2_WFC_2 :   Default_dir("wfpc2wf2.tab", table_name); break;
	   case WFPC2_WFC_3 :   Default_dir("wfpc2wf3.tab", table_name); break;
	   case WFPC2_WFC_4 :   Default_dir("wfpc2wf4.tab", table_name); break;
	       case FOC_F48 :   
	       case FOC_F96 :   Default_dir("foc.tab", table_name); break;
	      case FOC2_F48 :
	      case FOC2_F96 :   Default_dir("foc2ab.tab", table_name); break;
              case NICMOS_1 :   Default_dir("nicmos1.tab", table_name); break;
              case NICMOS_2 :   Default_dir("nicmos2.tab", table_name); break;
              case NICMOS_3 :   Default_dir("nicmos3.tab", table_name); break;
           case NEWNICMOS_1 :   if ( Pars.old_nicmos == 0 )
					Default_dir("newnicmos1.tab", table_name); 
				else
					Default_dir("newnicmos1_old.tab", table_name); 
				break;
           case NEWNICMOS_2 :   if ( Pars.old_nicmos == 0 )
					Default_dir("newnicmos2.tab", table_name); 
				else
					Default_dir("newnicmos2_old.tab", table_name); 
				break;
           case NEWNICMOS_3 :   Default_dir("newnicmos3.tab", table_name); break;
	      case ACS_WFC1 :   Default_dir("acswfc1.tab", table_name); break;
	      case ACS_WFC2 :   Default_dir("acswfc2.tab", table_name); break;
	       case ACS_HRC :   Default_dir("acshrc.tab", table_name); break;
       case ACS_HRC_OFFSPOT :   Default_dir("acshrc.tab", table_name); break;
	       case ACS_SBC :   Default_dir("acssbc.tab", table_name); break;
	      case STIS_CCD :   Default_dir("stisccd.tab", table_name); break;
	      case STIS_NUV :   Default_dir("stisnuv.tab", table_name); break;
	      case STIS_FUV :   Default_dir("stisfuv.tab", table_name); break;
	    case WFC3_UVIS1 :   Default_dir("wfc3_uvis1.tab", table_name); break;
	    case WFC3_UVIS2 :   Default_dir("wfc3_uvis2.tab", table_name); break;
	       case WFC3_IR :   Default_dir("wfc3_ir.tab", table_name); break;
	}

	strcpy( Pars.zernike_name, table_name );

	if ( (file = fopen(table_name, "r")) == NULL )
	{
		fprintf(stderr, "Get_aberrations: Could not open %s\n", table_name );
		exit(2);
	}


	if ( Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3 )
		Nicmos_aberrations( file, obs_date );
	else if ( Pars.chip == NEWNICMOS_1 )
	{
		Pars.best_focus_pam = 1.755;
		Pars.observation_pam = 1.755;
		Read_zernike_data( file );
	}
	else if ( Pars.chip == NEWNICMOS_2 )
	{
		Pars.best_focus_pam = 0.143;
		Pars.observation_pam = 0.143;
		Read_zernike_data( file );
	}
	else if ( Pars.chip == NEWNICMOS_3 )
	{
		Pars.best_focus_pam = -11.06;
		Pars.observation_pam = -9.49;
		Read_zernike_data( file );
	}
	else if ( Pars.chip >= WFPC2_PC ) 
	{
		Read_zernike_data( file );
	}	
	else
		Set_old_aberrations( obs_date, file );

	fclose( file );

} /* Get_aberrations */

