/*  File     :  polypsf.c
 *
 *  Contents :
 *      Add_mono_to_poly  : weighted addition of a monochromatic PSF to the
 *                              current polychromatic PSF array.
 *      Compute_poly_psf : Compute polychromatic PSF for a given position.
 *
 *  Author   :  John Krist (STScI)
 *  Date     :  January 1992
 *
 *  Modified :
 *	March 1994 - JEK
 *	  Added PSF position to call to ComputeOPD in Computepoly_psf.
 *	  Added ConvolveScatter.
 *	  Added call to ConvolveScatter in Computepoly_psf.
 *
 *	August 1994 - JEK
 *	  Modified routines to use complex structures.
 *
 *      October 1996 - Richard Hook
 *        Modified to handle "table" driven case where filters are
 *        not being used but there are several wavelengths.
 *
 *      December 1997 - JEK
 *	  Removed ConvolveScatter and replaced it with Convolve.  Added
 *	  Computeacs_kernel and modified Computepoly_psf to apply appropriate
 *	  scattering kernel to WFPC2, STIS CCD, or ACS CCD images.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include "tinytim.h"


/*-------------------------------------------------------------------------
*  Add_mono_to_poly :
*       weighted addition of the current monochromatic PSF to the poly_psf
*       array.
*
*  Inputs :
*       wave_index :  Current wavelength index
*       mono_psf   :  Monochromatic PSF image.
*       poly_psf   :  Polychromatic PSF image.
*     psf_x, psf_y :  Detector locations of PSF center (only used for
*		          ACS HRC, otherwise ignored)
*
*-------------------------------------------------------------------------*/
static void Add_mono_to_poly( int wave_index, float **mono_psf, float **poly_psf,
		int psf_x, int psf_y )
{
	int     x, y;
	float   weight;


	/* Flux normalize mono_psf, and then add it to poly_psf  *
	 * with weighting.				       */

	Flux_normalize( &mono_psf[0][0], Pars.int_dim, Pars.int_dim );

	/* If ACS HRC, add red halo (if wavelength > 600 nm) */

	if ( (Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT) && 
	     Pars.wavelength[wave_index] >= 0.6 && Pars.add_halo != 0 )
		Hrc_halo( mono_psf, Pars.wavelength[wave_index], psf_x, psf_y );

	if ( Pars.num_waves < 2 )      /* Only one wavelength */
	{

		for ( y = 0; y < Pars.int_dim; ++y )
			for ( x = 0; x < Pars.int_dim; ++x )
				poly_psf[y][x] = mono_psf[y][x];
	}
	else
	{
		weight = Pars.weight[wave_index];
		for ( y = 0; y < Pars.int_dim; ++y )
			for ( x = 0; x < Pars.int_dim; ++x )
				poly_psf[y][x] += weight * mono_psf[y][x];
	}

} /* Add_mono_to_poly */

/*-------------------------------------------------------------------------
*  Compute_poly_psf :
*       Compute the polychromatic PSF for the current position.
*
*  Inputs :
*       pos_index   :  Current position index
*       poly_psf    :  Polychromatic PSF image array.
*       mono_psf    :  Monochromatic PSF image array.
*
*  Outputs :
*       Returns the polychromatic PSF in poly_psf.
*
*  Modifications :
*	John Krist - March 1992 :
*	Shift pupil function to origin before computing FFT.  Does not
*	significantly affect the results, but is mathematically the
*	right thing to do.
*-------------------------------------------------------------------------*/
void Compute_poly_psf( int pos_index, float **poly_psf, float **mono_psf )
{
	int     i, crit_dim, previous_dim, x, y, stat, psf_x, psf_y;
	int     files_written, pupil_allocated, skip_psf;
	complex **pupil;
	float   **temp, max_weight, weight_limit=1e10;
	char	pupil_name[MAX_STRING];

	/* Define name for the temporary pupil/OPD file */

	strcpy( pupil_name, Pars.root_name );
	strcat( pupil_name, ".pup" );

	psf_x = Pars.x[pos_index];
	psf_y = Pars.y[pos_index];

	for ( y = 0; y < 3; ++y )
		for ( x = 0; x < 3; ++x )
			Pars.weighted_kernel[y][x] = 0.0;

	for ( y = 0; y < Pars.int_dim; ++y )
		for ( x = 0; x < Pars.int_dim; ++x )
			poly_psf[y][x] = 0.0;

	/* if using smart_skip, determine max weight */

	if ( Pars.smart_skip != 0 )
	{
		max_weight = Pars.weight[0];
		for ( i = 1; i < Pars.num_waves; ++i )
		    if ( Pars.weight[i] > max_weight )
			max_weight = Pars.weight[i];
		weight_limit = Pars.weight_limit * max_weight;
	}

	/* Compute monochromatic PSF for each wavelength */

	files_written = 0;
	previous_dim = 0;
	pupil_allocated = 0;
	pupil = NULL;

	for ( i = 0; i < Pars.num_waves; ++i )
	{
		/* If the current critical psf dimension is not the same *
		 * as the previous one, or this is the first wavelength  *
		 * to be computed, then compute the pupil and OPD for    *
		 * the new dimension.  If the new arrays can be used for *
		 * succeeding wavelengths, write them out to a file so   *
		 * that they will not need to be recomputed.             */

                if ( Pars.num_waves > 1 )
		{
		   if ( (Pars.smart_skip != 0) && (Pars.weight[i] < weight_limit) )
		   {
			printf("   Skipping PSF %d/%d (low weight), wavelength = %.2f nm\n",
				i+1, Pars.num_waves, Pars.wavelength[i] * 1000. );
			skip_psf = 1;
		   }
		   else
		   {
	    	  	printf( "   Computing PSF %d/%d for wavelength %.2f nm (weight=%f)\n", 
		    		i+1, Pars.num_waves, Pars.wavelength[i] * 1000.,
		    		Pars.weight[i] );
			skip_psf = 0;
		   }
	        }
                else
		{
                   printf("   Computing PSF for wavelength = %.2f nm\n",
        		Pars.wavelength[i] * 1000.);
		   skip_psf = 0;
		}

		if ( skip_psf == 0 )
		{
		  crit_dim = Pars.crit_dim[i];

		  if ( crit_dim != previous_dim )
		  {
			if ( pupil_allocated )
				Free_complex_image( pupil );

			pupil = Alloc_complex_image( crit_dim, crit_dim );
			pupil_allocated = 1;

			temp = Alloc_image( crit_dim, crit_dim );

			Compute_pupil( crit_dim, psf_x, psf_y, temp );

			for ( y = 0; y < crit_dim; ++y )
			{
				for ( x = 0; x < crit_dim; ++x )
				{
					pupil[y][x].r = temp[y][x];
					temp[y][x] = 0.0;
				}
			}

			Compute_opd( crit_dim, temp, psf_x, psf_y );

			for ( y = 0; y < crit_dim; ++y )
				for ( x = 0; x < crit_dim; ++x )
					pupil[y][x].i = temp[y][x];

			Free_image( temp );

			Shift_complex_to_origin( pupil, crit_dim );

			/* if next wavelength has same grid size as the *
			 * current one, then save the pupil/OPD arrays  *
			 * to a temporary file so they don't have to be *
			 * recomputed.					*/

			if ( Pars.crit_dim[i+1] == crit_dim )
			{
				if ( files_written )
					Delete_file( pupil_name );

				stat = Write_image( pupil_name, pupil[0],
					crit_dim, crit_dim, COMPLEX );

				if ( stat == ERROR )
				{
					fprintf(stderr, "Error writing temp file\n");
					exit(2);
				}
				files_written = 1;
			}
		  }
		  else
		  {
			stat = Read_image( pupil_name, pupil[0], crit_dim, crit_dim, 
				COMPLEX );
			if ( stat == ERROR )
			{
				fprintf(stderr, "Error reading temp file\n");
				exit(2);
			}
		  }

		  previous_dim = crit_dim;

		  Compute_mono_psf( i, mono_psf, pupil );

		  /* if STIS CCD, compute the pixel charge diffusion                 *
		   * kernel for the current wavelength.  If the PSF is normally      *
		   * sampled, convolve it with the kernel.  If it is sub-sampled,    *
		   * the kernel is not applied; however, a final kernel is generated *
		   * that is the sum of the wavelength-weighted kernels.  This       *
		   * kernel is printed to the .log file.  The user should apply it   *
		   * after rebinning the subsampled PSF to normal sampling.  (Yeah,  *
		   * I know it's not really mathematically correct, but it's better  *
		   * than nothing.)						     */

		  if ( Pars.chip == STIS_CCD )
		  {
	    		Compute_acs_kernel( Pars.wavelength[i], Pars.weight[i] );
	    		if ( Pars.scatter_flag != 0 && Pars.psf_scale == Pars.pixel_size )
		  		Convolve( mono_psf, Pars.kernel, Pars.int_dim, Pars.int_dim );
		  }
		
		  /* Add, with weighting, the monochromatic PSF to the *
		   * polychromatic PSF array.                          */

		  Add_mono_to_poly( i, mono_psf, poly_psf, psf_x, psf_y );
	       }
	}

	Free_complex_image( pupil );

	if ( files_written )
		Delete_file( pupil_name );

	/* Apply pixel scattering and F1042M halo, if WFPC2 */

	if ( Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4 )
	{
	    if ( Pars.scatter_flag != 0 && Pars.psf_scale == Pars.pixel_size )
	    {
		/* Pars.kernel already defined for WFPC2 */
	        Convolve( poly_psf, Pars.kernel, Pars.int_dim, Pars.int_dim );
		printf("   Applied WFPC2 pixel scatter function.\n");
	    }

	    if ( strncasecmp(Pars.filter_name, "F1042M", 5) == 0 && Pars.add_halo > 0 )
	    {
		printf("adding F1042M CCD diffuse halo...\n");
		Halo_f1042m( poly_psf, Pars.int_dim, Pars.sampling );
	    }
	}

} /* Compute_poly_psf */
