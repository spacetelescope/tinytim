/* File      : rdpupil.c
 *
 * Contents  :
 *      Read_ota_data      : Read OTA pupil data
 *	Read_wfpc_data     : Read WFPC pupil data
 *	Read_wfpc2_data    : Read WFPC2 pupil data
 *	Read_foc_data      : Read FOC pupil data
 *	Read_nicmos_data   : Read NICMOS pupil data
 *      Read_wfc3uvis_data : Read WFC3 UVIS data
 *      Read_wfc3ir_data   : Read WFC3 IR data
 *	Read_acs_data	   : Read ACS pupil data
 *	Read_stis_data     : Read STIS pupil data
 *      Open_pupil_table   : Open a pupil table file for a specified camera
 *	Read_pupil_data    : Read pupil data for a specified camera
 * 	Read_pupil_table   : Open and read a pupil table for a specified camera
 *      Read_geom_data     : Read geometric information (aperture position, and geometric distortion) for a specified camera
 *
 * Author    : John Krist (STScI)
 * Date      : May 1993
 *
 * Modifications :
 *
 *	March 1994 - JEK
 *	  Modified Read_wfpc2_data to read in x and y focus slopes.
 *
 *      March 1996 - RNH
 *        Added NICMOS
 *
 *      April 1996 - RNH
 *        Added NICMOS cold mask obscurations
 *
 *      June 1997 - JEK
 *	  Added field and focus dependent aberration stuff to NICMOS
 *
 *      December 1997 - JEK
 *	  Added ACS, STIS pupil table read routines
 *
 *      December 1998 - JEK
 *	  Added WFPC2 field dependent aberration coefficient stuff
 *
 *	October 1999 - JEK
 *	  Removed pixel size from NICMOS pupil file - now read from another file
 *    
 *      January 2003 - JEK
 *        Added new NICMOS+cryocooler
 *
 *      March 2008 - RNH & FS
 *        Added WFC3 support - and moved distortion reading to a separate routine
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "tinytim.h"

#define  OLD_NICMOS  0
#define  NEW_NICMOS  1

/*--------------------------------------------------------------------------*/
static void Read_ota_data( FILE *file )
{
	int	i;

	/* Read OTA secondary and spider dimensions */

	sscanf( Get_entry( file ), "%f", &Pars.ota_secondary_radius );
	sscanf( Get_entry( file ), "%f", &Pars.ota_spider_width );

	/* Read OTA mirror pad positions and radii */

	for ( i = 0; i <= 2; ++i )
		sscanf( Get_entry( file ), "%f %f %f", 
			&Pars.ota_pad_v3[i], &Pars.ota_pad_v2[i], 
			&Pars.ota_pad_radius[i] );

} /* Read_ota_data */

/*--------------------------------------------------------------------------*/
static void Read_wfpc_data( FILE *file, int camera )
{
	/* Read wavelength range */

	sscanf( Get_entry( file ), "%f %f",
		&Pars.min_wavelength, &Pars.max_wavelength );

	/* Read WFPC alignment center */

	sscanf( Get_entry( file ), "%f %f", 
		&Pars.obsc_center_x, &Pars.obsc_center_y );

	/* Read WFPC secondary and spider dimensions */

	sscanf( Get_entry( file ), "%f", &Pars.camera_secondary_radius );
	sscanf( Get_entry( file ), "%f", &Pars.camera_spider_width );

	/* Read WFPC camera rotation (relative to OTA) */

	sscanf( Get_entry( file ), "%f", &Pars.theta );

	/* Read pixel size in arcsec */

	sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

	/* Read WFC entrance baffle info */

	if ( camera >= 1 && camera <= 4 )
	{
		sscanf( Get_entry( file ), "%f", &Pars.camera_entrance_radius );
		sscanf( Get_entry( file ), "%f %f %f %f",
				&Pars.wfc_vig_radius, &Pars.wfc_vig_xc,
				&Pars.wfc_vig_yc, &Pars.wfc_vig_inside );
	}

	/* Read detector center v2,v3 pupil coordinates */

	sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3 );

} /* Read_wfpc_data */

/*---------------------------------------------------------------------------*/
static void Read_wfpc2_data( FILE *file )
{
	int	i, j;

	/* Read wavelength range */

	sscanf( Get_entry( file ), "%f %f",
		&Pars.min_wavelength, &Pars.max_wavelength );

	/* Read WFPC2 obscuration alignment center */

	sscanf( Get_entry( file ), "%f %f", 
		&Pars.obsc_center_x, &Pars.obsc_center_y );

	/* Read WFPC2 secondary and spider dimensions */

	sscanf( Get_entry( file ), "%f", &Pars.camera_secondary_radius );
	sscanf( Get_entry( file ), "%f", &Pars.camera_spider_width);

	/* Read WFPC2 camera rotation (relative to OTA) */

	sscanf( Get_entry( file ), "%f", &Pars.theta );

	/* Read WFPC2 pixel size in arcsec */

	sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

	/* Read detector center v2,v3 pupil coordinates */

	sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3 );

	/* Read WFPC2 field dependent aberration coefficients */

	for ( j = 0; j < 5; ++j )        /* aberration index */
		sscanf( Get_entry( file ), "%g %g %g %g %g", 
			&Pars.caber[j][0], 
			&Pars.caber[j][1], 
			&Pars.caber[j][2], 
			&Pars.caber[j][3], 
			&Pars.caber[j][4] );

	/* Read WFPC2 pixel scattering kernel */

	for ( i = 0; i < 3; ++i )
		sscanf( Get_entry( file ), "%f %f %f", 
			&Pars.kernel[i][0],
			&Pars.kernel[i][1],
			&Pars.kernel[i][2] );

	for ( j = 0; j < 3; ++j )
	   for ( i = 0; i < 3; ++i )
		Pars.weighted_kernel[j][i] = Pars.kernel[j][i];
	
} /* Read_wfpc2_data */

/*---------------------------------------------------------------------------*/
static void Read_foc_data( FILE *file )
{
	/* Read wavelength range */

	sscanf( Get_entry( file ), "%f %f",
		&Pars.min_wavelength, &Pars.max_wavelength );

	/* Read FOC camera rotation (relative to OTA) */

	sscanf( Get_entry( file ), "%f", &Pars.theta );

	/* Read FOC pixel size */

	sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

	/* Read FOC axial offset */

	sscanf( Get_entry( file ), "%f %f", &Pars.v2,
		&Pars.v3 );

} /* Read_foc_data */

/*---------------------------------------------------------------------------*/
static void Read_nicmos_data( FILE *file, int camera, int nicmos_type )
{
	int	j;

	/* Read wavelength range */

	sscanf( Get_entry( file ), "%f %f",
		&Pars.min_wavelength, &Pars.max_wavelength );

      	/* Read NICMOS obscuration sizes */

      	sscanf( Get_entry( file ), "%f", &Pars.camera_secondary_radius );
      	sscanf( Get_entry( file ), "%f", &Pars.camera_spider_width );
      	sscanf( Get_entry( file ), "%f", &Pars.camera_entrance_radius );

	if ( nicmos_type == NEW_NICMOS )
	{
      		sscanf( Get_entry( file ), "%f", &Pars.nicmos_mask_rotation );
      		sscanf( Get_entry( file ), "%f", &Pars.ota_offset );
      		sscanf( Get_entry( file ), "%f %f", &Pars.nicmos_mask_x_off, &Pars.nicmos_mask_y_off );
	}

        /* Read NICMOS camera rotation */

        sscanf( Get_entry( file ), "%f", &Pars.theta );

        /* Read NICMOS axial offset */

        sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3);

	if ( nicmos_type == NEW_NICMOS )
		sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

        /* If camera 2 or 3, read the NICMOS cold-mask pad cover  *
         * parameters.  Camera 1 doesn't have these.              */

	if ( camera != NICMOS_1 && camera != NEWNICMOS_1 )
	{
		/* read v2,v3 coordinates of pads */

        	sscanf( Get_entry( file ), "%f %f %f", 
                	&Pars.nic_pad_v2[0],
                	&Pars.nic_pad_v2[1],
                	&Pars.nic_pad_v2[2]);

        	sscanf( Get_entry( file ), "%f %f %f", 
                	&Pars.nic_pad_v3[0],
                	&Pars.nic_pad_v3[1],
                	&Pars.nic_pad_v3[2]);

		/* read x,y sizes of pads, as proportions of the pupil radius */

        	sscanf( Get_entry( file ), "%f %f %f", 
                	&Pars.nic_pad_x[0],
                	&Pars.nic_pad_x[1],
                	&Pars.nic_pad_x[2]);

        	sscanf( Get_entry( file ), "%f %f %f", 
                	&Pars.nic_pad_y[0],
                	&Pars.nic_pad_y[1],
                	&Pars.nic_pad_y[2]);

		/* read rotations of pads */

        	sscanf( Get_entry( file ), "%f %f %f", 
                	&Pars.nic_pad_rot[0],
                	&Pars.nic_pad_rot[1],
                	&Pars.nic_pad_rot[2]);
	}

	/* read scaling factor to convert defocus in PAM space (mm) to *
	 * wavefront error (microns)                                   */

	sscanf( Get_entry( file ), "%f", &Pars.pam_to_z4 );

	/* read focus-dependent aberration matrices (microns) */

	sscanf( Get_entry( file ), "%f %f", 
		&Pars.z5_vs_z4[0], &Pars.z5_vs_z4[1] );
	sscanf( Get_entry( file ), "%f %f", 
		&Pars.z6_vs_z4[0], &Pars.z6_vs_z4[1] );
	sscanf( Get_entry( file ), "%f %f", 
		&Pars.z7_vs_z4[0], &Pars.z7_vs_z4[1] );
	sscanf( Get_entry( file ), "%f %f", 
		&Pars.z8_vs_z4[0], &Pars.z8_vs_z4[1] );
	sscanf( Get_entry( file ), "%f %f", 
		&Pars.z11_vs_z4[0], &Pars.z11_vs_z4[1] );

	/* read the field-dependent aberration matrices (microns) */

	for ( j = 0; j <= 2; ++j )
	 	sscanf( Get_entry( file ), "%f %f %f",
			&Pars.z4_vs_xy[j][0], &Pars.z4_vs_xy[j][1],
			&Pars.z4_vs_xy[j][2] );
	for ( j = 0; j <= 2; ++j )
	 	sscanf( Get_entry( file ), "%f %f %f",
			&Pars.z5_vs_xy[j][0], &Pars.z5_vs_xy[j][1],
			&Pars.z5_vs_xy[j][2] );
	for ( j = 0; j <= 2; ++j )
	 	sscanf( Get_entry( file ), "%f %f %f",
			&Pars.z6_vs_xy[j][0], &Pars.z6_vs_xy[j][1],
			&Pars.z6_vs_xy[j][2] );
	for ( j = 0; j <= 2; ++j )
	 	sscanf( Get_entry( file ), "%f %f %f",
			&Pars.z7_vs_xy[j][0], &Pars.z7_vs_xy[j][1],
			&Pars.z7_vs_xy[j][2] );
	for ( j = 0; j <= 2; ++j )
	 	sscanf( Get_entry( file ), "%f %f %f",
			&Pars.z8_vs_xy[j][0], &Pars.z8_vs_xy[j][1],
			&Pars.z8_vs_xy[j][2] );

} /* Read_nicmos_data */

/*---------------------------------------------------------------------------*/
static void Read_wfc3ir_data( FILE *file )
{
        int i,j;

	/* Read WFC3 IR channel wavelength range */

	sscanf( Get_entry( file ), "%f %f",
		&Pars.min_wavelength, &Pars.max_wavelength );

      	/* Read WFC3 IR channel obscuration sizes */

      	sscanf( Get_entry( file ), "%f", &Pars.camera_secondary_radius );
      	sscanf( Get_entry( file ), "%f", &Pars.camera_spider_width );
      	sscanf( Get_entry( file ), "%f", &Pars.camera_entrance_radius );

	sscanf( Get_entry( file ), "%f %f", &Pars.wfc3_mask_x,
		&Pars.wfc3_mask_y );

	/* Read camera rotation */

        sscanf( Get_entry( file ), "%f", &Pars.theta );
	
        /* Read pixel size */

        sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

        /* Read camera axial offset */

        sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3);

        /* Read field dependent aberration coefficients.  Aberrations are, *
         * in order: focus, x & y astigmatism, x & y coma                  */

        for ( i = 0; i < 5; ++i )
            for ( j = 0; j < 6; ++j )
                sscanf( Get_entry( file ), "%f %f %f %f %f %f",
                        &Pars.wfc3_aber[i][j][0],
                        &Pars.wfc3_aber[i][j][1],
                        &Pars.wfc3_aber[i][j][2],
                        &Pars.wfc3_aber[i][j][3],
                        &Pars.wfc3_aber[i][j][4],
                        &Pars.wfc3_aber[i][j][5] );

        /* read geometric information - position of apertures and distortion coefficients */
        Read_geom_data( file);

        /* Read the charge diffusion kernels, and their spatial variation */
        for ( i = 0; i < NUM_WFC3_KERNELS; ++i )
        {
                sscanf( Get_entry( file ), "%f",
                        &Pars.wfc3_kernel_wavelength[i] );

                for ( j = 0; j < 3; ++j )
                        sscanf( Get_entry( file ), "%f %f %f",
                                &Pars.wfc3_kernel[i][j][0],
                                &Pars.wfc3_kernel[i][j][1],
                                &Pars.wfc3_kernel[i][j][2]);
        }

        sscanf( Get_entry(file), "%d", &Pars.blur_xy_num_wavelengths );

        for ( i = 0; i < Pars.blur_xy_num_wavelengths; ++i )
        {
                sscanf( Get_entry(file), "%f", &Pars.blur_xy_wavelength[i] );

                for ( j = 0; j <= 5; ++j )
                        sscanf( Get_entry(file), "%f %f %f %f %f %f",
                                &Pars.blur_xy_sigma[i][j][0], &Pars.blur_xy_sigma[i][j][1],
                                &Pars.blur_xy_sigma[i][j][2], &Pars.blur_xy_sigma[i][j][3],
                                &Pars.blur_xy_sigma[i][j][4], &Pars.blur_xy_sigma[i][j][5] );

        }

} /* Read_wfc3ir_data */

/*----------------------------------------------------------------------------*/
static void Read_wfc3uvis_data( FILE *file )
{
        int     i, j;

        /* Read wavelength range in nm */

        sscanf( Get_entry( file ), "%f %f",
                &Pars.min_wavelength, &Pars.max_wavelength );

        /* Read camera rotation (relative to OTA) */

        sscanf( Get_entry( file ), "%f", &Pars.theta );

        /* Read pixel size */

        sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

        /* Read field dependent aberration coefficients.  Aberrations are, *
         * in order: focus, x & y astigmatism, x & y coma                  */

        for ( i = 0; i < 5; ++i )
            for ( j = 0; j < 6; ++j )
                sscanf( Get_entry( file ), "%f %f %f %f %f %f",
                        &Pars.wfc3_aber[i][j][0],
                        &Pars.wfc3_aber[i][j][1],
                        &Pars.wfc3_aber[i][j][2],
                        &Pars.wfc3_aber[i][j][3],
                        &Pars.wfc3_aber[i][j][4],
                        &Pars.wfc3_aber[i][j][5] );

        /* Read geometric information - position of apertures and distortion coefficients */
        Read_geom_data( file );

        /* Read the charge diffusion kernels, and their spatial variation */
        for ( i = 0; i < NUM_WFC3_KERNELS; ++i )
        {
                sscanf( Get_entry( file ), "%f",
                        &Pars.wfc3_kernel_wavelength[i] );

                for ( j = 0; j < 3; ++j )
                        sscanf( Get_entry( file ), "%f %f %f",
                                &Pars.wfc3_kernel[i][j][0],
                                &Pars.wfc3_kernel[i][j][1],
                                &Pars.wfc3_kernel[i][j][2]);
        }

        sscanf( Get_entry(file), "%d", &Pars.blur_xy_num_wavelengths );

        for ( i = 0; i < Pars.blur_xy_num_wavelengths; ++i )
        {
                sscanf( Get_entry(file), "%f", &Pars.blur_xy_wavelength[i] );

                for ( j = 0; j <= 5; ++j )
                        sscanf( Get_entry(file), "%f %f %f %f %f %f",
                                &Pars.blur_xy_sigma[i][j][0], &Pars.blur_xy_sigma[i][j][1],
                                &Pars.blur_xy_sigma[i][j][2], &Pars.blur_xy_sigma[i][j][3],
                                &Pars.blur_xy_sigma[i][j][4], &Pars.blur_xy_sigma[i][j][5] );

        }


} /* Read_wfc3uvis_data */

/*----------------------------------------------------------------------------*/
static void Read_acs_data( FILE *file, int camera )
{
	int	i, j;

	/* Read wavelength range in nm */

	sscanf( Get_entry( file ), "%f %f", 
	  	&Pars.min_wavelength, &Pars.max_wavelength );

	/* Read ACS camera rotation (relative to OTA) */

	sscanf( Get_entry( file ), "%f", &Pars.theta );

	/* Read pixel size */

	sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

	if ( camera == ACS_HRC_OFFSPOT )
	{
		sscanf( Get_entry( file ), "%f", &Pars.acs_outer_radius );
		sscanf( Get_entry( file ), "%f", &Pars.acs_secondary_radius );
		sscanf( Get_entry( file ), "%f", &Pars.acs_spider_width );
		sscanf( Get_entry( file ), "%f", &Pars.acs_pad_radius );
		sscanf( Get_entry( file ), "%f %f", &Pars.acs_stop_shift, &Pars.acs_stop_angle );
	}

	/* Read field dependent aberration coefficients.  Aberrations are, *
	 * in order: focus, x & y astigmatism, x & y coma	           */

	for ( i = 0; i < 5; ++i )
	    for ( j = 0; j < 6; ++j )
		sscanf( Get_entry( file ), "%f %f %f %f %f %f",
			&Pars.acs_aber[i][j][0],
			&Pars.acs_aber[i][j][1],
			&Pars.acs_aber[i][j][2],
			&Pars.acs_aber[i][j][3],
			&Pars.acs_aber[i][j][4],
			&Pars.acs_aber[i][j][5] );

        /* Read the geometric information */
        Read_geom_data( file );

	/* Read CCD pixel scattering kernels */

	if ( camera == ACS_HRC || camera == ACS_HRC_OFFSPOT || camera == ACS_WFC1 || camera == ACS_WFC2 )
	{
	   	for ( i = 0; i < NUM_ACS_KERNELS; ++i )
	   	{
			sscanf( Get_entry( file ), "%f", 
				&Pars.acs_kernel_wavelength[i] );

			for ( j = 0; j < 3; ++j )
		    		sscanf( Get_entry( file ), "%f %f %f", 
					&Pars.acs_kernel[i][j][0],
					&Pars.acs_kernel[i][j][1],
					&Pars.acs_kernel[i][j][2]);
	   	}

		sscanf( Get_entry(file), "%d", &Pars.blur_xy_num_wavelengths );
		for ( i = 0; i < Pars.blur_xy_num_wavelengths; ++i )
		{
			sscanf( Get_entry(file), "%f", &Pars.blur_xy_wavelength[i] );
			for ( j = 0; j <= 5; ++j )
		    		sscanf( Get_entry(file), "%f %f %f %f %f %f",
					&Pars.blur_xy_sigma[i][j][0], &Pars.blur_xy_sigma[i][j][1], 
					&Pars.blur_xy_sigma[i][j][2], &Pars.blur_xy_sigma[i][j][3], 
					&Pars.blur_xy_sigma[i][j][4], &Pars.blur_xy_sigma[i][j][5] );
		}
	}

} /* Read_acs_data */

/*---------------------------------------------------------------------------*/
static void Read_stis_data( FILE *file, int camera )
{
	int	i, j;

	/* Read wavelength range */

	sscanf( Get_entry( file ), "%f %f", &Pars.min_wavelength, 
		&Pars.max_wavelength );

	/* Read STIS camera rotation (relative to OTA) */

	sscanf( Get_entry( file ), "%f", &Pars.theta );

	/* Read STIS pixel size */

	sscanf( Get_entry( file ), "%f", &Pars.pixel_size );

	/* Read STIS axial offsets */

	sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3 );

	if ( camera == STIS_CCD )
	{
		/* Read STIS Lyot stop parameters (CCD only) */

		sscanf( Get_entry( file ), "%f %f", 
			&Pars.stis_stop_x_radius, &Pars.stis_stop_y_radius );
		sscanf( Get_entry( file ), "%f %f", 
			&Pars.stis_stop_x_offset, &Pars.stis_stop_y_offset );

		/* Read CCD pixel scattering kernels (put into ACS kernel array) */

		for ( i = 0; i < NUM_ACS_KERNELS; ++i )
		{
			sscanf( Get_entry( file ), "%f", 
				&Pars.acs_kernel_wavelength[i] );

			for ( j = 0; j < 3; ++j )
		    		sscanf( Get_entry( file ), "%f %f %f", 
					&Pars.acs_kernel[i][j][0],
					&Pars.acs_kernel[i][j][1],
					&Pars.acs_kernel[i][j][2]);
		}
	}

} /* Read_stis_data */

/*---------------------------------------------------------------------------*/
static FILE *Open_pupil_table( int camera, char *table_name )
{
	FILE	*file;

	switch ( camera )
	{
		         case   1  : Default_dir( "wfpcwf1.pup", table_name );
				     break;
		         case   2  : Default_dir( "wfpcwf2.pup", table_name );
				     break;
		         case   3  : Default_dir( "wfpcwf3.pup", table_name );
				     break;
		         case   4  : Default_dir( "wfpcwf4.pup", table_name );
				     break;
		         case   5  : Default_dir( "wfpcpc5.pup", table_name );
				     break;
		         case   6  : Default_dir( "wfpcpc6.pup", table_name );
				     break;
		         case   7  : Default_dir( "wfpcpc7.pup", table_name );
				     break;
		         case   8  : Default_dir( "wfpcpc8.pup", table_name );
				     break;
    
		     case WFPC2_PC : Default_dir( "wfpc2pc1.pup", table_name );
				     break;
		  case WFPC2_WFC_2 : Default_dir( "wfpc2wf2.pup", table_name );
				     break; 
		  case WFPC2_WFC_3 : Default_dir( "wfpc2wf3.pup", table_name );
				     break; 
		  case WFPC2_WFC_4 : Default_dir( "wfpc2wf4.pup", table_name );
				     break; 

		     case FOC_F48  : Default_dir( "focf48.pup", table_name );
                                     break;
		     case FOC_F96  : Default_dir( "focf96.pup", table_name );
				     break;

		     case FOC2_F48 : Default_dir( "costarf48.pup", table_name );
				     break;
		     case FOC2_F96 : Default_dir( "costarf96.pup", table_name );
				     break;

                     case NICMOS_1 : Default_dir( "nicmos1.pup", table_name );
				     break;
                     case NICMOS_2 : Default_dir( "nicmos2.pup", table_name );
				     break;
                     case NICMOS_3 : Default_dir( "nicmos3.pup", table_name );
                                     break;

                  case NEWNICMOS_1 : Default_dir( "newnicmos1.pup", table_name );
				     break;
                  case NEWNICMOS_2 : if ( Pars.old_nicmos == 0 )
					Default_dir( "newnicmos2.pup", table_name );
				     else
					Default_dir( "newnicmos2_old.pup", table_name );
				     break;
                  case NEWNICMOS_3 : Default_dir( "newnicmos3.pup", table_name );
                                     break;

                     case ACS_WFC1 : Default_dir( "acswfc1.pup", table_name );
				     break;
                     case ACS_WFC2 : Default_dir( "acswfc2.pup", table_name );
				     break;
		      case ACS_HRC : Default_dir( "acshrc.pup", table_name );
				     break;
	      case ACS_HRC_OFFSPOT : Default_dir( "acshrcoffspot.pup", table_name );
				     break;
		      case ACS_SBC : Default_dir( "acssbc.pup", table_name );
				     break;

		     case STIS_CCD : Default_dir( "stisccd.pup", table_name );
                                     break;
		     case STIS_NUV : Default_dir( "stisnuv.pup", table_name );
                                     break;
		     case STIS_FUV : Default_dir( "stisfuv.pup", table_name );
				     break;

		   case WFC3_UVIS1 : Default_dir( "wfc3_uvis1.pup", table_name );
				     break;

		   case WFC3_UVIS2 : Default_dir( "wfc3_uvis2.pup", table_name );
				     break;

		     case WFC3_IR :  Default_dir( "wfc3_ir.pup", table_name );
				     break;

			   default : fprintf(stderr, "**Unknown Camera %d**\n",
					camera );
				     exit(2);
	}

	if ( (file = fopen( table_name, "r" )) == NULL )
	{
		fprintf(stderr, "Open_pupil_table : Could not open %s\n",
			table_name );
		exit(2);
	}

	return( file );

} /* Open_pupil_table */

/*----------------------------------------------------------------------------*/
void Read_pupil_data( int camera, FILE *file )
{

	Read_ota_data( file );

	if ( camera >= 1 && camera <= 8 )
		Read_wfpc_data( file, camera );
	else if ( camera >= WFPC2_PC && camera <= WFPC2_WFC_4 )
		Read_wfpc2_data( file );
	else
	{
		switch ( camera )
		{
		     case FOC2_F48 :
		     case FOC2_F96 : 
		     case FOC_F48  :
		     case FOC_F96  : Read_foc_data( file );
				     break;

                     case NICMOS_1 :
                     case NICMOS_2 :
                     case NICMOS_3 : Read_nicmos_data( file, camera, OLD_NICMOS );
                                     break;

                  case NEWNICMOS_1 :
                  case NEWNICMOS_2 :
                  case NEWNICMOS_3 : Read_nicmos_data( file, camera, NEW_NICMOS );
                                     break;

		     case ACS_WFC1 :
		     case ACS_WFC2 :
		     case  ACS_HRC :
	      case ACS_HRC_OFFSPOT :
		     case  ACS_SBC : Read_acs_data( file, camera );
				     break;

		     case STIS_CCD :
		     case STIS_NUV :
		     case STIS_FUV : Read_stis_data( file, camera );
				     break;

	           case WFC3_UVIS1 : 
                   case WFC3_UVIS2 : Read_wfc3uvis_data( file );
				     break;

		     case WFC3_IR  : Read_wfc3ir_data( file );
				     break;

			   default : fprintf(stderr, "**Unknown Camera %d**\n",
					camera );
				     exit(2);
		}
	}

} /* Read_pupil_data */

/*----------------------------------------------------------------------------*/
void Read_pupil_table( int camera, char *table_name )
{
	FILE	*file;

	file = Open_pupil_table( camera, table_name );
	Read_pupil_data( camera, file );
	fclose( file );

} /* Read_pupil_table */

/*----------------------------------------------------------------------------*/
void Read_geom_data( FILE *file)
{

       /* Read axial offsets for reference position */

        sscanf( Get_entry( file ), "%f %f", &Pars.v2, &Pars.v3 );
        
        /* Read detector pixel reference positions */
        
        sscanf( Get_entry( file ), "%f %f", &Pars.x_ref, &Pars.y_ref );
                        
        /* Read V2,V3 pupil coordinates of camera (not detector) center */

        sscanf( Get_entry( file ), "%f %f", &Pars.v2c, &Pars.v3c );

        /* Read field X,Y to v2,v3 transform coefficients;  NOTE: these *
         * coefficients are relative to the camera optical axis         */

        sscanf( Get_entry(file), "%f %f", &Pars.xy_to_xc[0], &Pars.xy_to_xc[1] );
        sscanf( Get_entry(file), "%f %f %f",
                &Pars.xy_to_xc[2], &Pars.xy_to_xc[3], &Pars.xy_to_xc[4] );
        sscanf( Get_entry(file), "%f %f %f %f", 
                &Pars.xy_to_xc[5], &Pars.xy_to_xc[6], &Pars.xy_to_xc[7], &Pars.xy_to_xc[8] );
        sscanf( Get_entry(file), "%f %f %f %f %f",
                &Pars.xy_to_xc[9], &Pars.xy_to_xc[10], &Pars.xy_to_xc[11], &Pars.xy_to_xc[12], &Pars.xy_to_xc[13] );
        
        sscanf( Get_entry(file), "%f %f", &Pars.xy_to_yc[0], &Pars.xy_to_yc[1] );
        sscanf( Get_entry(file), "%f %f %f",
                &Pars.xy_to_yc[2], &Pars.xy_to_yc[3], &Pars.xy_to_yc[4] );
        sscanf( Get_entry(file), "%f %f %f %f",
                &Pars.xy_to_yc[5], &Pars.xy_to_yc[6], &Pars.xy_to_yc[7], &Pars.xy_to_yc[8] );
        sscanf( Get_entry(file), "%f %f %f %f %f",
                &Pars.xy_to_yc[9], &Pars.xy_to_yc[10], &Pars.xy_to_yc[11], &Pars.xy_to_yc[12], &Pars.xy_to_yc[13] );

        /* Read v2,v3 to field X,Y transform coefficients */

        sscanf( Get_entry(file), "%f %f", &Pars.xcyc_to_x[0], &Pars.xcyc_to_x[1] );
        sscanf( Get_entry(file), "%f %f %f",
                &Pars.xcyc_to_x[2], &Pars.xcyc_to_x[3], &Pars.xcyc_to_x[4] );
        sscanf( Get_entry(file), "%f %f %f %f",
                &Pars.xcyc_to_x[5], &Pars.xcyc_to_x[6], &Pars.xcyc_to_x[7], &Pars.xcyc_to_x[8] );
        sscanf( Get_entry(file), "%f %f %f %f %f",
                &Pars.xcyc_to_x[9], &Pars.xcyc_to_x[10], &Pars.xcyc_to_x[11], &Pars.xcyc_to_x[12], &Pars.xcyc_to_x[13] );

        sscanf( Get_entry(file), "%f %f", &Pars.xcyc_to_y[0], &Pars.xcyc_to_y[1] );
        sscanf( Get_entry(file), "%f %f %f",
                &Pars.xcyc_to_y[2], &Pars.xcyc_to_y[3], &Pars.xcyc_to_y[4] );
        sscanf( Get_entry(file), "%f %f %f %f",
                &Pars.xcyc_to_y[5], &Pars.xcyc_to_y[6], &Pars.xcyc_to_y[7], &Pars.xcyc_to_y[8] );
        sscanf( Get_entry(file), "%f %f %f %f %f",
                &Pars.xcyc_to_y[9], &Pars.xcyc_to_y[10], &Pars.xcyc_to_y[11], &Pars.xcyc_to_y[12], &Pars.xcyc_to_y[13] );

} /* Read_geom_data */
