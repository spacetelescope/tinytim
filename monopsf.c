/* File      :  monopsf.c
 *
 * Contents  :
 *      WaveCorrectOPD :  Convert OPD from 547 nm to specified wavelength.
 *      Computemono_psf :  Compute the monochromatic PSF.
 *
 * Author    :  John Krist (STScI)
 * Date      :  January 1992
 *
 * Modifications :
 *	JEK - May 1992
 *		Added code to use NewIntegratePsf.
 *
 *      JEK - December 1992
 *		Added code to write out the pupil and wavefront maps if
 *		the Writepupil or WriteWave flags are set.
 *
 *	JEK - August 1994
 *		Modified routines to use complex structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tinytim.h"

#define APERTURE   1
#define WAVEFRONT  2

/*-------------------------------------------------------------------------
*  WaveCorrectOPD :
*       Convert OPD from wave at 0.547 microns to waves at wavelength microns.
*
*  Inputs :
*       opd_dim     : dimension of pupil function array.
*       wavelength : wavelength in microns of desired OPD.
*       pupil      : Complex valued pupil function defined at 547 nm.
*
*  Returns :
*       Returns the wavelength corrected OPD in the pupil array.
*-----------------------------------------------------------------------------*/
static void Wave_correct_opd( int opd_dim, float wavelength, complex **pupil )
{
	int     x, y;
	float   c;

	c = 0.547 / wavelength;

	for ( y = 0; y < opd_dim; ++y )
		for ( x = 0; x < opd_dim; ++x )
			pupil[y][x].i *= c;

} /* Wave_correct_opd */

/*-----------------------------------------------------------------------------
*  Write_pupil_map :
*       Write the aperture or opd file to disk.
*
*  Inputs :
*	image    :  image array to be written to disk.
*	dim      :  Size of image array (dim by dim).
*	name     :  String containing name of array, either "pupil" or "wave"
*       data_type :  Flag indicating pupil (1) or opd (2) data.
*
*----------------------------------------------------------------------------*/
static void Write_pupil_map( complex **image, int dim, char *name, int data_type )
{
	float	**temp;
	int	i, j;
	char    temp_name[MAX_STRING];


	/* Create a temporary array for the shifted image */

	temp = Alloc_image( dim, dim );

	if ( data_type == APERTURE )
	{
		for ( j = 0; j < dim; ++j )
			for ( i = 0; i < dim; ++i )
				temp[j][i] = image[j][i].r;
	}
	else
	{
		for ( j = 0; j < dim; ++j )
			for ( i = 0; i < dim; ++i )
				temp[j][i] = image[j][i].i;
	}

	Shift_to_center( temp, dim );

       	printf("  ** Writing %s map to tinytim_%s.fits ...", name, name );
        fflush( stdout );
	sprintf( temp_name, "tinytim_%s.fits", name );
        Write_FITS( temp_name, temp[0], dim, dim, MAP_FILE, 0 );

	printf(" Done\n");
	fflush( stdout );

	Free_image( temp );

} /* Write_pupil_map */

/*-----------------------------------------------------------------------------
*  Compute_mono_psf :
*       Compute the monochromatic PSF for a specified wavelength.
*
*  Inputs :
*       w       :  wavelength index
*      mono_psf :  Array to contain monochromatic PSF.
*       pupil   :  pupil function array (complex structure).
*
*  Returns :
*       Returns the monochromatic, detector integrated PSF in mono_psf.
*----------------------------------------------------------------------------*/
void Compute_mono_psf( int w, float **mono_psf, complex **pupil )
{
	int	size_in, size_out, i, j;
	float	**temp, scale_in, scale_out, wavelength;


	size_in = Pars.crit_dim[w];

	/* Convert OPD to current wavelength */

	Wave_correct_opd( size_in, Pars.wavelength[w], pupil );

	/* If the Writepupil or WriteWave flags are set, write the arrays *
	 * to disk.							  */

	if ( Pars.write_pupil == 1 )
	{
		Write_pupil_map( pupil, size_in, "pupil", APERTURE );
		Pars.write_pupil = 0;
	}

	if ( Pars.write_wave == 1 )
	{
		Write_pupil_map( pupil, size_in, "wave", WAVEFRONT );
		Pars.write_wave = 0;
	}

	Compute_crit_psf( pupil, size_in );   /* pupil now has critical PSF */

	/* If Pars.write_crit_psf = 1, then write the critically sampled  *
	 * PSF to the file critpsf.fits.  The PSF is first shifted to    *
	 * the center.  The program then stops (hoping that the exit()   *
	 * routine will take care of cleaning up the arrays).		 */

	if ( Pars.write_crit_psf == 1 )
	{
		Shift_complex_to_center( pupil, size_in );

		temp = Alloc_image( size_in, size_in );

		for ( j = 0; j < size_in; ++j )
			for ( i = 0; i < size_in; ++i )
				temp[j][i] = pupil[j][i].r;

	        printf("*** Writing critical PSF to critpsf.fits ***\n");
	        Write_FITS( "critpsf.fits", temp[0], size_in, size_in, PSF_FILE, 1 );

		Free_image( temp );

		printf("Stopping.\n");
		exit(0);
	}

	/* Integrate PSF */

        wavelength = Pars.wavelength[w];
        scale_in = 360. * 3600. * wavelength / (4 * M_PI * DIAM_HST_MICRONS);

        size_out = Pars.int_dim;
        scale_out = Pars.psf_scale;

	Integrate_psf( pupil, scale_in, size_in, mono_psf, scale_out, size_out, 
		Pars.jitter_major, Pars.jitter_minor, Pars.jitter_angle );

} /* Compute_mono_psf */

