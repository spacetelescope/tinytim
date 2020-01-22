/* File       :  param.c
 *
 * Contents   :
 *	Read_param_data : used by Read_params()
 *	Read_params : Read parameter file.
 *
 * Author     :  John Krist (STScI)
 * Date       :  January 1992
 *
 * Modified   :
 *	 Aug 1992 - J.K.
 *		Removed input for Wfpc2PCRot and replaced NumSincSamples
 *		with SincLimit (which is a floating point number).
 *
 *	 Dec 1992 - J.K.
 *		Added lines to read in pupil/wavemap flags.
 *
 *	 March 1994 - J.K.
 *		Added lines to read scattering and focus slope flags.
 *
 *       June 1997 - J.K.
 *              Added NICMOS focus and XY dependent aberration flags 
 *
 *	 Dec 1997 - J.K.
 *		Added STIS, ACS.
 *
 *	 Oct 1998 - J.K.
 *		Modified filter parameters.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tinytim.h"

/*----------------------------------------------------------------------*/
static void Read_param_data( FILE *file )
{
	char	temp[MAX_STRING];
	int	i, x, y;

	for ( y = 0; y < 3; ++y )
		for ( x = 0; x < 3; ++x )
			Pars.weighted_kernel[y][x] = 0.0;

	/* First line must begin with &BeginParn.n */

	sprintf( temp, "&BeginPar%3.1f", (float)VERSION_NUM );

	if ( strstr(Get_entry( file), temp) == NULL )
	{
		fprintf(stderr, "\nERROR! Pre-Version %3.1f parameter file?\n",
			(float)VERSION_NUM);
		exit(2);
	}

	sscanf( Get_entry(file), "%s", Pars.root_name );
	sscanf( Get_entry(file), "%s", Pars.zernike_name );
	sscanf( Get_entry(file), "%s", Pars.pupil_name );
	Get_entry(file);  /* Skip date */
	sscanf( Get_entry(file), "%s", Pars.camera_name );
	sscanf( Get_entry(file), "%d", &Pars.chip );
	sscanf( Get_entry(file), "%s", Pars.filter_name );
	sscanf( Get_entry(file), "%s", Pars.spectrum_file );
	sscanf( Get_entry(file), "%f", &Pars.jitter_major );
	sscanf( Get_entry(file), "%f", &Pars.jitter_minor );
	sscanf( Get_entry(file), "%f", &Pars.jitter_angle );

	sscanf( Get_entry(file), "%d", &Pars.num_pos );
	for ( i = 0; i < Pars.num_pos; ++i )
		sscanf( Get_entry(file), "%d %d", 
			&Pars.x[i], &Pars.y[i] );

	sscanf( Get_entry(file), "%f", &Pars.psf_size );

	sscanf( Get_entry(file), "%d", &Pars.smart_skip );
	sscanf( Get_entry(file), "%f", &Pars.weight_limit );

	sscanf( Get_entry(file), "%d", &Pars.num_waves );

	for ( i = 0; i < Pars.num_waves; ++i )
		sscanf( Get_entry(file), "%f %f %d", 
			&Pars.wavelength[i], &Pars.weight[i], &Pars.crit_dim[i] );

	sscanf( Get_entry(file), "%d", &Pars.int_dim );
	sscanf( Get_entry(file), "%f", &Pars.psf_scale );
	sscanf( Get_entry(file), "%f", &Pars.sampling );
	sscanf( Get_entry(file), "%d", &Pars.use_map );

	Pars.adjust_for_xy = 0;
	if ( (Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4) || 
	     (Pars.chip >= ACS_WFC1 && Pars.chip <= ACS_SBC) ||
             (Pars.chip >= WFC3_UVIS1 && Pars.chip <= WFC3_IR) )
	{
		sscanf( Get_entry(file), "%d", &Pars.adjust_for_xy );
	}

	Pars.scatter_flag = 0;	
	if ( (Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4) || 
	     Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 || Pars.chip == ACS_HRC ||
	     Pars.chip == ACS_HRC_OFFSPOT || Pars.chip == STIS_CCD ||
             Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 || Pars.chip == WFC3_IR )
	{
		sscanf( Get_entry(file), "%d", &Pars.scatter_flag );
	}

	Pars.add_halo = 0;
	if ( (Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4) ||
             (Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT) )
	{
		sscanf( Get_entry(file), "%d", &Pars.add_halo );
	}

	if ( Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3 )
	{
		sscanf( Get_entry(file), "%f", &Pars.best_focus_pam );
		sscanf( Get_entry(file), "%f", &Pars.observation_pam );
		sscanf( Get_entry(file), "%d", &Pars.adjust_for_focus );
		sscanf( Get_entry(file), "%d", &Pars.adjust_for_xy );
		sscanf( Get_entry(file), "%f", &Pars.ota_offset );
		sscanf( Get_entry(file), "%f", &Pars.nicmos_mask_rotation );
		sscanf( Get_entry(file), "%f %f", &Pars.nicmos_mask_x_off, &Pars.nicmos_mask_y_off );
	}

	if ( Pars.chip >= NEWNICMOS_1 && Pars.chip <= NEWNICMOS_3 )
	{
		sscanf( Get_entry(file), "%f", &Pars.best_focus_pam );
		sscanf( Get_entry(file), "%f", &Pars.observation_pam );
		sscanf( Get_entry(file), "%d", &Pars.adjust_for_focus );
		sscanf( Get_entry(file), "%d", &Pars.adjust_for_xy );
	}

	sscanf( Get_entry(file), "%d", &Pars.write_pupil );
	sscanf( Get_entry(file), "%d", &Pars.write_wave );
	sscanf( Get_entry(file), "%d", &Pars.write_crit_psf );

} /* Read_param_data */

/*----------------------------------------------------------------------
*  Read_params :
*	Read input parameters from a file.
*
*  Inputs :
*	param_file  :  Name of parameter file.
*
*  Returns :
*	Fills parameter structure.
*----------------------------------------------------------------------*/
void Read_params( char *param_file )
{
	FILE	*file;

	if ( (file = fopen( param_file, "r" )) == NULL )
	{
		fprintf(stderr, "Error opening %s\n", param_file );
		exit(2);
	}

	Init_entry_list();   /* Needed because Get_entry writes to list */

	Read_param_data( file );
	Read_pupil_data( Pars.chip, file );
	Read_zernike_data( file );

	Delete_entry_list();

	fclose( file );

} /* Read_params */

