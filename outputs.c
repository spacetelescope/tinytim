/* File      : outputs.c
 *
 * Contents  :
 *      Write_params       : Write Tiny Tim parameters to a file.
 *
 * Author    : John Krist (STSci)
 * Date      : January 1992
 *
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "tinytim.h"

/*----------------------------------------------------------------------
*  Write_params :
*	Write output parameters to a file.
*
*  Inputs :
*	param_file :  Name of parameter file.
*       obs_date   :  Pointer to a DateStruct observation date.
*----------------------------------------------------------------------*/
void Write_params( char *param_file, struct DateStruct *obs_date )
{
	FILE	*file;
	int	i;

	if ( (file = fopen( param_file, "w" )) == NULL )
	{
		fprintf(stderr, "Error opening %s\n", param_file );
		exit(2);
	}

	fprintf( file, "&BeginPar%3.1f\n", (float)VERSION_NUM );

	fprintf( file, "%s  # Output PSF file rootname\n", Pars.root_name );
	fprintf( file, "%s  # Name of Zernike table file\n", Pars.zernike_name );
	fprintf( file, "%s  # Name of pupil table file\n", Pars.pupil_name );
	fprintf( file, "%.2d-%.2d-%.4d # Observation date (dd-mm-yyyy)\n",
		obs_date->day, obs_date->month, obs_date->year );
	fprintf( file, "%s  # Camera\n", Pars.camera_name );
	fprintf( file, "%d  # Camera ID number\n", Pars.chip );
	fprintf( file, "%s  # Filter\n", Pars.filter_name );
	fprintf( file, "%s  # Spectrum file\n", Pars.spectrum_file );
	fprintf( file, "%f  # Major axis jitter in mas (-1 = no jitter)\n", Pars.jitter_major);
	fprintf( file, "%f  # Minor axis jitter in mas\n", Pars.jitter_minor);
	fprintf( file, "%f  # Angle of jitter major axis in deg from +X axis\n", Pars.jitter_angle);
	fprintf( file, "%d  # Number of positions\n", Pars.num_pos );
	for ( i = 0; i < Pars.num_pos; ++i )
		fprintf( file, "%d %d  # Position %d\n", Pars.x[i], Pars.y[i], i+1 );

	fprintf( file, "%.3f  # PSF diameter in arcsecs\n", Pars.psf_size );

	fprintf( file, "%d  # Skip wavelengths with low weights? (0=no)\n",
		Pars.smart_skip );
	fprintf( file, "%f  # Min good weight limit coefficient\n",
		Pars.weight_limit );

	fprintf( file, "%d  # Number of wavelengths\n", Pars.num_waves );

	for ( i = 0; i < Pars.num_waves; ++i )
		fprintf( file, "%f %f %d # Wavelength %d (microns), weight, grid size\n", 
			Pars.wavelength[i], Pars.weight[i], Pars.crit_dim[i], i+1 );

	fprintf( file, "%d  # Integrated PSF dimension (pixels)\n", Pars.int_dim );
	fprintf( file, "%f  # Integrated PSF scaling (arcsec)\n", Pars.psf_scale );
	fprintf( file, "%f  # Subsampling factor (1 = normal)\n", Pars.sampling );
	fprintf( file, "%d  #  Use mirror maps? (0 = no, otherwise yes)\n", 
	    Pars.use_map );

	if ( Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4 )
	{
		fprintf( file, "%d  #  Adjust for WFPC2 field aberrations? (0=no)\n",
			Pars.adjust_for_xy );
		fprintf( file, "%d  #  Apply WFPC2 pixel scattering (0=no)\n",
			Pars.scatter_flag );
		fprintf( file, "1  #  Apply WFPC2 halo (if F1042M) (0=no)\n");
	}

	if ( Pars.chip >= ACS_WFC1 && Pars.chip <= ACS_SBC )
	{
		fprintf( file, "%d  #  Adjust for ACS field aberrations? (0=no)\n",
			Pars.adjust_for_xy );
		if ( Pars.chip != ACS_SBC )
			fprintf( file, "%d  #  Apply ACS pixel scattering (0=no)\n",
				Pars.scatter_flag );
		if ( Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT )
			fprintf( file, "1  #  Apply HRC red halo (0=no)\n" );
	}

        if ( Pars.chip >= WFC3_UVIS1 && Pars.chip <= WFC3_IR )
        {
                fprintf( file, "%d  #  Adjust for WFC3 field aberrations? (0=no)\n",
                        Pars.adjust_for_xy );
                fprintf( file, "%d  #  Apply WFC3 pixel scattering (0=no)\n",
                        Pars.scatter_flag );
        }

	if ( Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3 )
	{
		fprintf( file, "%f # PAM position of camera best focus (mm)\n",
			Pars.best_focus_pam );
		fprintf( file, "%f # PAM position of observation (mm)\n",
			Pars.observation_pam );
		fprintf( file, "%d # Adjust aberrations w.r.t. focus? (0=no)\n",
			Pars.adjust_for_focus );
		fprintf( file, "%d # Adjust aberrations w.r.t. X,Y? (0=no)\n",
			Pars.adjust_for_xy );
		fprintf( file, "%f # OTA offset rotation\n", Pars.ota_offset );
		fprintf( file, "%f # NICMOS cold mask rotation\n",
			Pars.nicmos_mask_rotation );
		fprintf( file, "%f %f # NICMOS cold mask V2,V3 offsets\n",
			Pars.nicmos_mask_x_off, Pars.nicmos_mask_y_off );
	}

	if ( Pars.chip >= NEWNICMOS_1 && Pars.chip <= NEWNICMOS_3 )
	{
		fprintf( file, "%f # PAM position of camera best focus (mm)\n",
			Pars.best_focus_pam );
		fprintf( file, "%f # PAM position of observation (mm)\n",
			Pars.observation_pam );
		fprintf( file, "%d # Adjust aberrations w.r.t. focus? (0=no)\n",
			Pars.adjust_for_focus );
		fprintf( file, "%d # Adjust aberrations w.r.t. X,Y? (0=no)\n",
			Pars.adjust_for_xy );
	}

	if ( Pars.chip == STIS_CCD )
	{
		fprintf( file, "%d  #  Apply STIS pixel scattering (0=no)\n",
			Pars.scatter_flag );
	}

	fprintf( file, "0  # Write out pupil map at 1st wavelength?\n" );
	fprintf( file, "0  # Write out wave map at 1st wavelength?\n" );
	fprintf( file, "0  # Write out crit psf at 1st wavelength & stop?\n" );

	/* WriteEntryList writes out pupil and Zernike parameters */

	Write_entry_list( file );

	fclose( file );

} /* Write_params */

