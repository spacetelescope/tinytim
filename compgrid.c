/*  File     :  compgrid.c (formerly inputs3.c)
 *
 *  Contents :
 *	Compute_grid_dims : Compute the PSF grid dimensions
 *
 *  Author   :  John Krist (STScI)
 *  Date     :  January 1992
 *
 *  Modifications :
 *	J. Krist - March 1992
 *	   Changed ComputeGridDims to use the focal ratios and pixel sizes
 *	   from the pupil table.
 *
 *      J. Krist - June 1993
 *	   Got rid of maximum sampling factor, since users can now choose
 *	   whatever sampling they like.  Changed algorithm to choose the
 *	   critically sampled PSF grid sizes so that the critically sampled
 *	   PSF is at least 11 arcsec for WFPC2 or COSTAR and 8 sec for the
 *	   others.
 *
 *	J. Krist - July 1994
 *	   Modified ComputeGridDims to allow 15" or 30" diameter WFPC1 & WFPC2
 *	   PSFs if the wavelength is greater than 300 nm or 550 nm, 
 *	   respectively.
 *
 *      R. N. Hook - March 1996
 *         Modified to include NICMOS.
 *
 *      J. Krist - June 1997
 *	   Fixed bug in NICMOS code (broke on large PSFs).  Added printout
 *         of assumed f/ratio and pixel size.
 *
 *	J. Krist - Sept 1997
 *	   Set maximum PSF size for NICMOS to 30 or less, rather than just 30
 *
 *	J. Krist - Dec 1997
 *	   Added STIS & ACS; cleaned up code a little
 */

#include <stdio.h>
#include <math.h>
#include "tinytim.h"

/* Size of array containing allowable critical dimensions */
#define NUM_CRIT_DIMS  3

/*------------------------------------------------------------------------
*  Compute_grid_dims :
*	Compute the critical and integrated PSF dimensions.
*
*  Inputs:
*	chip	     :  chip / camera id.
*	num_waves    :  Number of wavelengths at (dimension of waves).
*	waves	     :  Array containing wavelengths at which to generate
*			  monochromatic PSFs (in microns).
*  Outputs:
*	psf_diameter :  Diameter of PSF in arcseconds.
*	crit_dim     :  Array variable that will contain the Nyquist-sampled PSF
*			  grid size (in pixels) for each wavelength
*	int_dim	     :  Variable which will contain the integrated PSF
*			  dimension (in pixels)
*-------------------------------------------------------------------------*/
void Compute_grid_dims( float *psf_diameter, int chip, int num_waves, float *waves, 
	int *crit_dim, int *int_dim )
{
	static  int crit_dim_table[NUM_CRIT_DIMS] = {512,1024,1280};
	int	i, j, max_grid_dim, max_int_dim;
	float	max_int_diameter, min_crit_diam, req_grid_dim, nyquist0, pixel_size;
	float	crit_scale;

	pixel_size = Pars.pixel_size;  /* nominal pixel size in arcsec */

	/* nyquist0 = nyquist sampling in arcsec at shortest wavelength */

	nyquist0 = 360.0 * 3600.0 * waves[0] / (4 * M_PI * DIAM_HST_MICRONS);

	/* Determine the largest integrated PSF size based on the shortest  *
	 * given wavelength and the largest possible critical grid size.    * 
	 * The maximum diameter for the FOC and COSTAR+FOC is 7.0".  For    *
	 * For WF/PC-1,WFPC2,ACS,& STIS CCD it's 15.0", if the wavelength   *
	 * is greater than 300 nm, or 30" if it's greater than 550 nm.  The *
	 * PSFs, however, are not accurate out that far. 		    */

	max_grid_dim = crit_dim_table[NUM_CRIT_DIMS-1];
	max_int_dim =  max_grid_dim * nyquist0 / pixel_size;
	max_int_dim = max_int_dim - max_int_dim % 2;      /* Make dim even */
	max_int_diameter = max_int_dim * pixel_size;
	if ( max_int_diameter > 30.0 )
		max_int_diameter = 30.0;
	max_int_diameter = floor(max_int_diameter*10.) / 10.0;

	if ( (chip < ACS_WFC1) || (chip > WFC3_IR))
	  printf("\nAssuming detector pixel size of %7.5f arcsec\n", pixel_size );
	else
	  printf("\nUsing undistorted critical sampling pixel size of %7.5f arcsec\n", pixel_size);

	do {
		printf("\nThe maximum computable PSF size is %.1f arcsec.\n",
			max_int_diameter);
		if ( chip >= WFPC2_PC )
		     printf("The recommended size is about 3.0 arcseconds.\n");
		else
		     printf("The recommended size is about 6.0 arcseconds.\n");
		printf("What diameter should your PSF be (in arcseconds)? : ");
		scanf("%f", psf_diameter); 
		*int_dim = (int)(*psf_diameter / pixel_size);
	} while ( *psf_diameter > max_int_diameter );

	if ( chip >= ACS_WFC1 && chip <= ACS_SBC )
		*int_dim += 4;

	/* Determine minimum critically sampled PSF grid diameter.  WFPC2  *
	 * and COSTAR+FOC PSFs should be computed using a critically       *
	 * sampled PSF of at least 11.0" diameter, while the others should *
	 * use at least 8.0".						   */

	if ( (chip <= 8) || (chip >= NICMOS_1 && chip <= NICMOS_3) 
	   || (chip >= NEWNICMOS_1 && chip <= NEWNICMOS_3) )
	{
		if ( *psf_diameter > 8.0 )
			min_crit_diam = *psf_diameter;
		else
			min_crit_diam = 8.0;
	}		
	else
	{
		if ( *psf_diameter > 11.0 )
			min_crit_diam = *psf_diameter;
		else
			min_crit_diam = 11.0;
	}

	/* Make the integrated dimension an even number (integration  *
	 * routine assumes it is even, otherwise the PSF peak would   *
	 * at the intersection of pixels.			      */

	*int_dim = *int_dim - *int_dim % 2;
	*psf_diameter = *int_dim * pixel_size;

	/* Search for grid size which will produce a critically sampled  *
	 * PSF which is large enough (>= MinCritDiam).  Determine this   *
	 * gridsize for each wavelength.				 */

	for ( i = 0; i < num_waves; ++i )
	{
		/* critscale = Nyquist sampling interval in arcsec */

		crit_scale = 360. * 3600. * waves[i] / (4 * M_PI * DIAM_HST_MICRONS);
		req_grid_dim = min_crit_diam / crit_scale;

		if ( req_grid_dim > crit_dim_table[NUM_CRIT_DIMS-1] )
		{
			crit_dim[i] = crit_dim_table[NUM_CRIT_DIMS-1];
		}
		else
		{
			for ( j = 0; j < NUM_CRIT_DIMS; ++j )
			{
				if ( req_grid_dim <= crit_dim_table[j] )
				{
					crit_dim[i] = crit_dim_table[j];
					break;
				}
			}
		}
	}

} /* Compute_grid_dims */

