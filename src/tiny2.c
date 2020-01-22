/*  File    :  tiny2.c (formerly mkpsf.c)
 *
 *  Contents: 
 *	main() : Main program for PSF generation routines.
 *
 *  Author  :  John Krist (STScI)
 *  Date    :  January 1992
 *
 *  Some modifications for WFC3
 *
 * Richard Hook & Felix Stoehr, ST-ECF, March 2008
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#ifndef MACOSX
#include <malloc.h>
#else
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "tinytim.h"

#ifdef TT_THREADED
#include <unistd.h>
#endif

static void TT3_param_file( char *tt3_param_file )
{
	FILE	*file;
	int	x_range, y_range;


	if ( Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 )
	{
		x_range = 4095;
		y_range = 2047;
	}
	else
	{
		x_range = 1023;
		y_range = 1023;
	}

        if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 )
        {
                x_range = 4095;
                y_range = 2050;
        }

        if ( Pars.chip == WFC3_IR )
        {  
                x_range = 1013;
                y_range = 1013;
        }

	file = fopen( tt3_param_file, "w" );

	fprintf( file, "# %s : optional parameter file for tiny3\n", tt3_param_file );
	fprintf( file, "#\n" );
	fprintf( file, "# Note : you do not have to modify this file to produce an ACS or WFC3 PSF\n" );
	fprintf( file, "# unless you want to generate a simulated observation convolved\n" );
	fprintf( file, "# with a PSF and then mapped onto a distorted grid.\n" );
	fprintf( file, "###################################################################\n" );
	fprintf( file, "INPUT_FILENAME=no_file \n" );
	fprintf( file, "PIXEL_SCALE=%f\n", Pars.psf_scale );
	fprintf( file, "X_CENTER=%d \n", Pars.x[0] );
	fprintf( file, "Y_CENTER=%d \n", Pars.y[0] );
	fprintf( file, "OUTPUT_FILENAME=no_file \n" ); 
	fprintf( file, "###################################################################\n" ); 
	fprintf( file, "# Required entries :\n" );
	fprintf( file, "#      INPUT_FILENAME\n" );
	fprintf( file, "#      PIXEL_SCALE\n" );
	fprintf( file, "#      X_CENTER (range is 0-%d)\n", x_range);
	fprintf( file, "#      Y_CENTER (range is 0-%d)\n", y_range);
	fprintf( file, "#      OUTPUT_FILENAME\n" );
	fprintf( file, "# Optional entries :\n" );
	fprintf( file, "#      PSF\n" );
	fprintf( file, "# Comment lines begin with '#'\n" );

	fclose( file );
}

/*---------------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
	int     pos, i, is_position_dependent, is_wfpc2, is_nicmos, is_acs, is_wfc3;
	float   **poly_psf, **mono_psf;
	char    fname[MAX_STRING], tt3_param_file[MAX_STRING];
	time_t  start_time, end_time;

#ifdef TT_THREADED
#ifdef MACOSX        
	int	mib[2], numproc;
	size_t  len;
#endif
#endif

	if ( argc < 2 )
	{
		fprintf(stderr, "Call is : tiny2 paramfile\n");
		exit(2);
	}


	printf("Tiny Tim v%3.1f", (float)VERSION_NUM);

	Num_threads = 1;

#ifdef TT_THREADED
#ifndef MACOSX
	Num_threads = sysconf(_SC_NPROCESSORS_ONLN);
#else
	mib[0] = CTL_HW;
	mib[1] = HW_NCPU;
	len = sizeof(int);
	sysctl(mib, 2, &numproc, &len, NULL, 0);
	Num_threads = numproc;
#endif
#endif
	if ( Num_threads > 1 )
   {
		printf(" (Multithreaded, %d CPUs detected)", Num_threads);
      if (Num_threads > MAX_THREADS)
      {
         printf(" limiting use to only %d CPUs at once",MAX_THREADS);
         Num_threads = MAX_THREADS;
      }
   }
	printf("\n");
 
	Read_params( argv[1] );  /* Read parameter file */

	is_wfpc2 = 0;
	is_nicmos = 0;
	is_acs = 0;
        is_wfc3 = 0;

	if ( Pars.chip >= WFPC2_PC && Pars.chip <= WFPC2_WFC_4 )
		is_wfpc2 = 1;
	else if ( Pars.chip >= NICMOS_1 && Pars.chip <= NICMOS_3 )
		is_nicmos = 1;
	else if ( Pars.chip >= ACS_WFC1 && Pars.chip <= ACS_SBC )
		is_acs = 1;
        else if ( Pars.chip >= WFC3_UVIS1 && Pars.chip <= WFC3_IR )
                is_wfc3 = 1;

	is_position_dependent = 0;
	if ( Pars.chip <= 8 || is_wfpc2 || is_nicmos || is_acs || Pars.chip == STIS_CCD || is_wfc3 )
		is_position_dependent = 1;

	/* Allocate image array for polychromatic and monochromatic PSFs */

	poly_psf = Alloc_image( Pars.int_dim, Pars.int_dim );
	mono_psf = Alloc_image( Pars.int_dim, Pars.int_dim );

	/* if the PSF is subsampled, then it will not be convolved with *
	 * the charge diffusion kernel (if WFPC2 or STIS CCD).  Put the *
	 * kernel in the FITS header so that the user can apply it      *
	 * after he has rebinned the PSF to normal sampling.  Ignore    *
	 * the ACS here as it is handled in tiny3.                      */

	Pars.num_comments = 0;

	time( &start_time );

	if ( is_acs || is_wfc3 )
		printf("Intermediate PSF dimensions are ");
	else
		printf("Integrated PSF dimensions are ");
	printf("%d by %d\n", Pars.int_dim, Pars.int_dim );

	if ( Pars.psf_scale != Pars.pixel_size && (is_wfpc2 || Pars.chip == STIS_CCD) ) 
	{
		Pars.num_comments = 6;
		for ( i = 0; i < Pars.num_comments; ++i )
			Pars.comment[i] = (char *)malloc(80);
	}

	/* Compute polychromatic PSF for each position */

	for ( pos = 0; pos < Pars.num_pos; ++pos )
	{
		Pars.current_pos = pos;
		if ( is_position_dependent )
			printf("\nComputing PSF for position %d/%d (x,y) = %d %d\n",
				pos+1, Pars.num_pos, Pars.x[pos], Pars.y[pos] );

		Compute_poly_psf( pos, poly_psf, mono_psf );

		if ( Pars.num_pos <= 100 )
		{
			if ( !is_acs && !is_wfc3 )
				sprintf( fname, "%s%02d.fits", Pars.root_name, pos );
			else
				sprintf( fname, "%s%02d_psf.fits", Pars.root_name, pos );
		}
		else
		{
			if ( !is_acs && !is_wfc3 )
				sprintf( fname, "%s%03d.fits", Pars.root_name, pos );
			else
				sprintf( fname, "%s%03d_psf.fits", Pars.root_name, pos );
		}

		printf("   Writing PSF to %s\n", fname );

		if ( Pars.psf_scale != Pars.pixel_size && (is_wfpc2 || Pars.chip == STIS_CCD) ) 
		{
			sprintf( Pars.comment[0], 
				"This PSF is subsampled, so the charge diffusion kernel was" ); 
			sprintf( Pars.comment[1],
				"NOT applied.  After you rebin this PSF to normal sampling," );
			sprintf( Pars.comment[2],
				"you should convolve it with the following kernel :" ); 
			for ( i = 0; i <= 2; ++i )
				sprintf( Pars.comment[i+3], "   %f %f %f\n", 
				 	Pars.weighted_kernel[i][0], 
				 	Pars.weighted_kernel[i][1], 
				 	Pars.weighted_kernel[i][2] );
		}

		Write_FITS( fname, poly_psf[0], Pars.int_dim, Pars.int_dim, PSF_FILE, 1 );
	}

	Free_image( poly_psf );
	Free_image( mono_psf );

	time( &end_time );
	printf("\nStarted at  %s", ctime( &start_time ));
	printf("Finished at %s", ctime( &end_time ));

	if ( Pars.num_comments != 0 )
		for ( i = 0; i < Pars.num_comments; ++i )
			free( Pars.comment[i] );

	/* if ACS or WFC3, write out template optional parameter file for tiny3 */

	if ( is_acs || is_wfc3)
	{
		sprintf( tt3_param_file, "%s.tt3", Pars.root_name );
		printf( "\nWriting template optional parameter file for tiny3 to %s.\n",
			tt3_param_file );
		TT3_param_file( tt3_param_file );

		printf("\nTo continue PSF processing for ACS and WFC3, you must run tiny3 to resample\n");
		printf("and distort the PSF.  You may also process a simulated scene (see\n");
		printf("the manual for details).\n\n" );
		printf("Just to distort the PSF, issue this command :\n");
		printf("\n        tiny3 %s\n\n", argv[1]); 
	}

	return(0);

} /* main */

