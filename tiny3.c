/*  File    :  tiny3.c 
 *
 *  Contents: 
 *	Tiny3_command_line : process command line arguments
 *	Read_psf : read PSF FITS file produced by tiny2
 *
 *  Author  :  John Krist
 *  Date    :  September 2000
 *
 *  Updated: Richard Hook & Felix Stoehr, March 2008
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

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


/*-------------------------------------------------------------------------*/
static void Tiny3_command_line( int argc, char *argv[], int *position, int *sub_factor, float *undistorted_scale )
{
	int	i;

	if ( argc < 2 )
	{
		fprintf(stderr, "Call is : tiny3 paramfile [POS=n] [SUB=n] [UNDIS=pixsize] \n");
		exit(0);
	}
	else if ( argc == 2 )
		return;

	for ( i = 2; i < argc; ++i )
	{
		if ( Find_string(argv[i], "POS=") != NULL )
		{
			sscanf( &argv[i][4], "%d", position );
			if ( *position < 0 || *position >= Pars.num_pos )
			{
				printf("ERROR : Position number must be 0 to %d\n", Pars.num_pos-1);
				exit(0);
			}
		} 
		else if ( Find_string(argv[i], "SUB=") != NULL )
		{
			sscanf( &argv[i][4], "%d", sub_factor );
			if ( *sub_factor < 1 || *sub_factor > 10 )
			{
				printf("ERROR : Subsampling factor must be an integer between 1-10\n");
				exit(0);
			}
		}
		else if ( Find_string(argv[i], "UNDIS=") != NULL )
		{
			sscanf( &argv[i][6], "%f", undistorted_scale );
			if ( *undistorted_scale < 0.0 )
			{
				printf("ERROR : Undistorted pixel scale must be >0 arcsec\n");
				exit(0);
			}
		}
		else if ( Find_string(argv[i], "=") != NULL )
		{
			printf( "ERROR : Unknown keyword (%s).\n", argv[i] );
			exit(0);
		}
		else
			strcpy( Pars.tt3_param_file, argv[i] ); 
	}
}

/*----------------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
	time_t  start_time, end_time;
	int	i, dx, dy, position, scene_nx, scene_ny, scene_nx_out, scene_ny_out;
	float	**psf_in, **psf_out, **scene_in, **convolved_scene, **scene_out, dist, min_dist=1e10;
	float	cscale, sigma, x_scale, y_scale;
	char	name[MAX_STRING];
                
#ifdef TT_THREADED
#ifdef MACOSX
        int	mib[2], numproc;
	size_t  len;
#endif
#endif

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
 
	position = 0;
	Pars.sub_factor = 1;

	Read_params( argv[1] );  /* Read parameter file */

	if ( Pars.chip < ACS_WFC1 || Pars.chip > WFC3_IR )
	{
		printf( "*** ERROR : Tiny3 currently only works on ACS & WFC3 PSFs ***\n" );
		exit(0);
	}

	Pars.tt3_param_file[0] = '\0';
	Pars.undistorted_scale = 0.0;
	Tiny3_command_line( argc, argv, &position, &Pars.sub_factor, &Pars.undistorted_scale );

	strcpy( Pars.scene_input_file, "no_file" );
	if ( Pars.tt3_param_file[0] != '\0' )
		Read_optional_parameters();


	if ( Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT )
	{
		x_scale = 0.028;
		y_scale = 0.025;
		cscale = 1.20;
	}
	else if ( Pars.chip == ACS_SBC )
	{
		x_scale = 0.0336;
		y_scale = 0.0301;
		cscale = 1.20;
	}
        else if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 )
        {
                x_scale = 0.040;
                y_scale = 0.040;
                cscale = 1.20;
        }
        else if ( Pars.chip == WFC3_IR )
        {
                x_scale = 0.13;
                y_scale = 0.13;
                cscale = 1.20;
        }
	else if (Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 )
	{
		x_scale = 0.0495;
		y_scale = 0.0495;
		cscale = 1.20;
	}
        else 
        {
                printf( "*** ERROR : Unknown chip/camera, tiny3 only works for ACS & WFC3 PSFs ***\n" );
                exit(0);
        }


	time( &start_time );

	Init_damped_interp();

        /* Compute the scattering convolution kernel */
        if ( Pars.scatter_flag ) 
        {
               if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 || Pars.chip == WFC3_IR ) 
                        for ( i = 0; i < Pars.num_waves; ++i )
          		        Compute_wfc3_kernel( Pars.wavelength[i], Pars.weight[i] );
               else 
                        for ( i = 0; i < Pars.num_waves; ++i )
                                Compute_acs_kernel( Pars.wavelength[i], Pars.weight[i] );
         }

 	if ( strcmp(Pars.scene_input_file, "no_file") == 0 )
 	{
		/* process PSF only; there is no input scene to deal with */

		Pars.pixel_size = Pars.pixel_size / Pars.sub_factor;

		Pars.distorted_dim = Pars.sub_factor * Pars.int_dim * Pars.psf_scale / y_scale * cscale;
		psf_out = Alloc_image( Pars.distorted_dim, Pars.distorted_dim );

		Pars.current_pos = position;
		printf("Processing PSF for position %d/%d : (x,y) = %d %d\n",
			position+1, Pars.num_pos, Pars.x[position], Pars.y[position] );

 		if ( Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 || Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT || 
                     Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 || Pars.chip == WFC3_IR )
			sigma = Compute_kernel_xy( Pars.x[position], Pars.y[position], Pars.wavelength[0] );

		psf_in = Read_psf( position );

		printf( "  Mapping PSF onto distorted grid.\n" );
		Distort( psf_in, Pars.int_dim, Pars.int_dim, Pars.x[position], Pars.y[position], 
		    psf_out, Pars.distorted_dim, Pars.distorted_dim, Pars.chip );

		Pars.num_comments = 0;

		if ( Pars.sub_factor == 1 && Pars.chip != ACS_SBC )
		{
			if ( Pars.scatter_flag )
			{
				printf( "  Convolving PSF with charge diffusion kernel.\n" );
				Convolve( psf_out, Pars.weighted_kernel, Pars.distorted_dim, Pars.distorted_dim );
				Convolve( psf_out, Pars.kernel_xy, Pars.distorted_dim, Pars.distorted_dim );
			}
			else
			{
				printf( "  NOTE : PSF not convolved with charge diffusion kernel.\n" );
			}
		} 
	        else if ( Pars.sub_factor != 1 && Pars.chip != ACS_SBC )
		{
			printf( "  NOTE : Subsampled, so not convolving with charge diffusion kernel.\n" );
                	Pars.num_comments = 6;
                	for ( i = 0; i < Pars.num_comments; ++i )
                        	Pars.comment[i] = (char *)malloc(80);

                	sprintf( Pars.comment[0],
                        	"This PSF is subsampled, so the charge diffusion kernel was" );
                	sprintf( Pars.comment[1],
                        	"NOT applied.  After you rebin this PSF to normal sampling," );
                	sprintf( Pars.comment[2],
                        	"you should convolve it with the following kernel :" );

			Convolve_kernel( Pars.weighted_kernel, Pars.kernel_xy );

                	for ( i = 0; i <= 2; ++i )
                        	sprintf( Pars.comment[i+3], "   %f %f %f",
                                	Pars.weighted_kernel[i][0],
                                 	Pars.weighted_kernel[i][1],
                                 	Pars.weighted_kernel[i][2] );
		}

		Flux_normalize( psf_out[0], Pars.distorted_dim, Pars.distorted_dim ); 

		/* compute aberrations for inclusion in FITS header */

		Compute_aberrations( Pars.x[position], Pars.y[position] );

		if ( Pars.num_pos <= 100 )
			sprintf( name, "%s%02d.fits", Pars.root_name, position );
		else
			sprintf( name, "%s%03d.fits", Pars.root_name, position );

		Pars.psf_scale = Pars.pixel_size;

		printf("  Writing distorted PSF to %s (%d by %d pixels)\n", 
			name, Pars.distorted_dim, Pars.distorted_dim );
		Write_FITS( name, psf_out[0], Pars.distorted_dim, Pars.distorted_dim, DISTORTED_PSF_FILE, 1 );

		Free_image( psf_in );
		Free_image( psf_out );

		if ( Pars.num_comments != 0 )
			for ( i = 0; i < Pars.num_comments; ++i )
				free( Pars.comment[i] );
	}
	else
	{
		/* an input scene has been specified, so PSF convolve and distort it */

		strcpy( name, Pars.scene_output_file );

		/* look for closest PSF position to scene center */

		for ( i = 0; i < Pars.num_pos; ++i )
		{
			dx = Pars.scene_x_center - Pars.x[i];
			dy = Pars.scene_y_center - Pars.y[i];
			dist = sqrt(dx*dx + dy*dy);
			if ( i == 0 )
				min_dist = dist;

			if ( dist <= min_dist )
			{
				min_dist = dist;
				position = i;
			}
		}

		if ( Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 || Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT ||
                     Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 || Pars.chip == WFC3_IR )
       
		{
			sigma = Compute_kernel_xy( Pars.x[position], Pars.y[position], Pars.wavelength[0] );
			Convolve_kernel( Pars.weighted_kernel, Pars.kernel_xy );
		}

		scene_in = Read_scene( &scene_nx, &scene_ny );

		/* if scene_in is NULL, it means that a blank image is used */

		if ( scene_in != NULL )
			psf_in = Read_psf( position );
		else
			psf_in = NULL;

		convolved_scene = Convolve_scene( scene_in, &scene_nx, &scene_ny, psf_in );

		if ( psf_in != NULL )
		{
			Free_image( psf_in );
			Free_image( scene_in );
		}

		if ( Pars.num_scene_psfs > 0 )
		{
			printf( "  Adding PSFs. " );
			psf_in = Read_psf( position );
			Add_psfs( psf_in, convolved_scene, scene_nx, scene_ny );
			Free_image( psf_in );
		}

		scene_nx_out = scene_nx * Pars.psf_scale / x_scale * cscale + 2;
		scene_ny_out = scene_ny * Pars.psf_scale / y_scale * cscale + 2;
		scene_out = Alloc_image( scene_nx_out, scene_ny_out );

		printf( "  Distorting scene (this may take a while).\n" );
		Distort( convolved_scene, scene_nx, scene_ny, Pars.scene_x_center, Pars.scene_y_center, 
			 scene_out, scene_nx_out, scene_ny_out, Pars.chip );

		if ( Pars.chip != ACS_SBC )
		{ 
			if ( Pars.scatter_flag == 1 )
			{
				printf( "  Convolving scene with charge diffusion kernel.\n" );
				Convolve( scene_out, Pars.weighted_kernel, scene_nx_out, scene_ny_out );
				Convolve( scene_out, Pars.kernel_xy, scene_nx_out, scene_ny_out );
			}
			else
			{
				printf( "  NOTE : scene not convolved with charge diffusion kernel.\n" );
			}
		}

		/* compute aberrations for inclusion in FITS header */

		Compute_aberrations( Pars.x[position], Pars.y[position] );

		printf( "  Writing distorted scene to %s (%d by %d pixels)\n", 
			name, scene_nx_out, scene_ny_out );
		Write_FITS( name, scene_out[0], scene_nx_out, scene_ny_out, IMAGE_FILE, 1 );

		Free_image( convolved_scene );
		Free_image( scene_out );
	}
 
	time( &end_time );
	printf("\nStarted at  %s", ctime( &start_time ));
	printf("Finished at %s", ctime( &end_time ));

	return(0);

} /* main */

