/* File : acs.c
 *
 * Contents :
 *	Read_psf : Read PSF produced by tiny2 from FITS file
 *	Read_scene : Read user-generated scene to be PSF convolved & distorted
 *	Convolve_scene : Resample, PSF convolve, and distort input scene 
 *	Resample_thread : Resample procedure used by Resample
 *	Resample : Resample input scene onto uniform grid
 *	Copy_to_center : Copy one image into the middle of another
 *	Read_optional_parameters : Read ACS input scene parameters from a file
 *	Add_psfs : Add PSFs to scene
 *
 * Written by John Krist, November 2000
 * Modified for WFC3 support, Richard Hook & Felix Stoehr, March 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tinytim.h"

#ifdef TT_THREADED
#include <pthread.h>
#endif

struct isparstruct {
	float	**image;
	int	nx, ny, dim_out;
	float	mag;
	complex **newimage;
};

struct sparstruct {
	struct isparstruct *ipars;
	int	thread;
	int	start_row, end_row;
};

/*--------------------------------------------------------------------------
 * Compute_kernel_xy :
 *    Compute the field-dependent charge diffusion kernel for the CCDs (ACS/WFC, WFC3/UVIS etc).
 *    Now also uses the same mechanism for WFC3/IR.
 *
 * Note - was called Compute_acs_kernel_xy.
 *
 * Renamed and generalised for WFC3, Richard Hook, March 2008.
 * WFC/IR added, Richard Hook, May 2008.
 *
 *------------------------------------------------------------------------*/
float Compute_kernel_xy( int xdet, int ydet, float lambda )
{
        int     i, iw, j;
        float   d, r, sigma, sigma0, sigma1, tot, w0, w1, x, y;


	if ( Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT || Pars.chip == WFC3_IR )
		d = 10.24;
	else
		d = 40.96;

        /* This offset is really only valid for ACS where an over-sized grid was used */
        x = xdet / d + 4;

	/* coefficients are for global x,y detector coordinates (assuming *
	 * that the chips are perfectly butted together)		  */

	if ( Pars.chip == ACS_WFC1 || Pars.chip == WFC3_UVIS1 )
		ydet = ydet + 2048;

        y = ydet / d + 4;

	/* compute sigma of gaussian kernel at wavelength "lambda" */

	for ( i = 0; i < Pars.blur_xy_num_wavelengths-1; ++i )
		if ( lambda < Pars.blur_xy_wavelength[i] )
			break;

	if ( i < Pars.blur_xy_num_wavelengths-1 )
	{
		iw = i;
		w0 = Pars.blur_xy_wavelength[iw];
		w1 = Pars.blur_xy_wavelength[iw+1];
                
        	sigma0 = 0.0;
        	sigma1 = 0.0;
        	for ( i = 0; i <= 5; ++i )
		{
                	for ( j = 0; j <= 5; ++j )
			{
                        	sigma0 += Pars.blur_xy_sigma[iw][i][j] * pow((float)x,(float)i) * pow((float)y,(float)j);
                        	sigma1 += Pars.blur_xy_sigma[iw+1][i][j] * pow((float)x,(float)i) * pow((float)y,(float)j);
                        }
		}
		sigma = (sigma1 - sigma0) / (w1 - w0) * (lambda - w0) + sigma0;
	}
	else
	{
		sigma = 0.0;
		iw = Pars.blur_xy_num_wavelengths - 1;
		for ( i = 0; i <= 5; ++i )
			for ( j = 0; j <= 5; ++j )
				sigma += Pars.blur_xy_sigma[iw][i][j] * pow((float)x,(float)i) * pow((float)y,(float)j);
	}

	/* Compute gaussian kernel */
        /* Added the case of sigma = 0.0 - just return a delta function */

        if ( sigma > 0.0 )
        {

      	      tot = 0.0;

	      for ( j = 0; j <= 2; ++j )
	      {
		     for ( i = 0; i <= 2; ++i )
		     {
		   	   r = sqrt((i-1)*(i-1)+(j-1)*(j-1));
			   Pars.kernel_xy[j][i] = exp(-0.5*(r/sigma)*(r/sigma));
			   tot += Pars.kernel_xy[j][i];
		     }
	       }

	       for ( j = 0; j <= 2; ++j )
	   	    for ( i = 0; i <= 2; ++i )
			Pars.kernel_xy[j][i] /= tot;
          }
          else
          {
             for ( j = 0; j <= 2; ++j )
                     for ( i = 0; i <= 2; ++i )
                           Pars.kernel_xy[j][i] = 0.0;

              Pars.kernel_xy[1][1] = 1.0;
           }

	return( sigma );

} /* Compute_kernel_xy */

/*-------------------------------------------------------------------------
*  Compute_acs_kernel :
*       Compute ACS WFC or HRC CCD pixel scatter (charge diffusion)
*       function at given wavelength.  Also used for STIS CCD.
*
*  Inputs :
*       wavelength : wavelength in microns
*       weight : PSF weight for current wavelength
*
*  Places kernel for current wavelength in Pars.kernel.  Keeps running weighted
*  sum of kernels in Pars.weighted_kernel.
*
*  Scatter function is obtained from tests on SITe CCDs with pinhole grids.
*  The function is known to have a subpixel dependence, but there is some
*  uncertainty in the measurements at that level.  The function varies with
*  wavelength, so three kernels are provided here (@ 0.4, 0.55, & 0.8 microns).
*  They are interpolated to derive the proper kernel for the input wavelength.
*  The kernels are 3 x 3 pixels at normal CCD pixel resolution and represent
*  the scattering function for an object placed at the center of a pixel.
*----------------------------------------------------------------------------*/
void Compute_acs_kernel( float wavelength, float weight )
{
        int 	i, iw, j;
	float	tot;

        iw = 0;
        while ( wavelength > Pars.acs_kernel_wavelength[iw] && iw < NUM_ACS_KERNELS-2 )
                ++iw;

	tot = 0.0;

        for ( j = 0; j <= 2; ++j )
        {
            for ( i = 0; i <= 2; ++i )
            {
                   Pars.kernel[j][i] =
                     (Pars.acs_kernel[iw+1][j][i] - Pars.acs_kernel[iw][j][i]) /
                     (Pars.acs_kernel_wavelength[iw+1] - Pars.acs_kernel_wavelength[iw])
                     * (wavelength - Pars.acs_kernel_wavelength[iw])
                     + Pars.acs_kernel[iw][j][i];
		   tot += Pars.kernel[j][i];
            }
        }

	for ( j = 0; j <= 2; ++j )
	    for ( i = 0; i <= 2; ++i )
		Pars.kernel[j][i] /= tot;

        for ( j = 0; j <= 2; ++j )
            for ( i = 0; i <= 2; ++i )
                   Pars.weighted_kernel[j][i] += (Pars.kernel[j][i] * weight);

} /* Compute_acs_kernel */

/*-------------------------------------------------------------------------
*  Compute_wfc3_kernel :
*       Compute WFC3 UVIS CCD pixel scatter (charge diffusion)
*       function at given wavelength.  Essentially the same as the ACS
*       version but separated to allow different logic in future.
*
*  Also used for the WFC3/IR channel - for the IPC effect rather than charge diffusion.
*
*  Note - the preliminary WFC2 UVIS charge diffusion kernels were at only
*         two wavelengths (250nm and 810nm), but the logic here work in this case too.
*
*   Richard Hook & Felix Stoehr, ST-ECF, March 2008.
*
*--------------------------------------------------------------------------*/

void Compute_wfc3_kernel( float wavelength, float weight )
{
        int     i, iw, j;
        float   tot;

        iw = 0;
        while ( wavelength > Pars.wfc3_kernel_wavelength[iw] && iw < NUM_WFC3_KERNELS-2 )
                ++iw;

        tot = 0.0;

        for ( j = 0; j <= 2; ++j )
        {
            for ( i = 0; i <= 2; ++i )
            {
                   Pars.kernel[j][i] =
                     (Pars.wfc3_kernel[iw+1][j][i] - Pars.wfc3_kernel[iw][j][i]) /
                     (Pars.wfc3_kernel_wavelength[iw+1] - Pars.wfc3_kernel_wavelength[iw])
                     * (wavelength - Pars.wfc3_kernel_wavelength[iw])
                     + Pars.wfc3_kernel[iw][j][i];

                   tot += Pars.kernel[j][i];
            }
        }

        for ( j = 0; j <= 2; ++j )
            for ( i = 0; i <= 2; ++i )
                Pars.kernel[j][i] /= tot;

        for ( j = 0; j <= 2; ++j )
            for ( i = 0; i <= 2; ++i )
                   Pars.weighted_kernel[j][i] += (Pars.kernel[j][i] * weight);

} /* Compute_wfc3_kernel */

/*--------------------------------------------------------------------------
* Convolve_kernel :
*  Convolve one 3 x 3 kernel with another, then replace "kernel1" with the result
*--------------------------------------------------------------------------*/
void Convolve_kernel( float kernel1[3][3], float kernel2[3][3] )
{
	float	k1[5][5], k2[5][5], tot;
	int	i, j, ii, jj;


	for ( j = 0; j <= 4; ++j )
	{
		for ( i = 0; i <= 4; ++i )
		{
			k1[j][i] = 0.0;
			k2[j][i] = 0.0;
		}
	}

	for ( j = 0; j <= 2; ++j )
		for ( i = 0; i <= 2; ++i )
			k1[j+1][i+1] = kernel1[j][i];

	for ( j = 0; j <= 4; ++j )
	{
		for ( i = 0; i <= 4; ++i )
		{
			for ( jj = -1; jj <= 1; ++jj )
			{
				if ( j+jj < 0 || j+jj > 4 )
					continue;

				for ( ii = -1; ii <= 1; ++ii )
				{
					if ( i+ii < 0 || i+ii > 4 )
						continue;

					k2[j][i] += k1[jj+j][ii+i] * kernel2[jj+1][ii+1];
				}
			}
		}
	}

	tot = 0.0;

	for ( j = 0; j <= 2; ++j )
	{
		for ( i = 0; i <= 2; ++i )
		{
			kernel1[j][i] = k2[j+1][i+1];
			tot += kernel1[j][i];
		}
	}

	for ( j = 0; j <= 2; ++j )
		for ( i = 0; i <= 2; ++i )
			kernel1[j][i] /= tot;

} /* Convolve_kernel */
	
/*--------------------------------------------------------------------------
* Read_psf :
*	Read the ACS or WFC3 subsampled PSF produced by tiny2 from a FITS file.
* Inputs :
*	position : Position index of PSF (0 is first)
* Outputs :
*	Returns image.  Note that this routine allocates the memory for
*	the image with Alloc_image().
*--------------------------------------------------------------------------*/
float **Read_psf( int position )
{
        int     nx, ny;
        float   **psf;

	if ( Pars.num_pos <= 100 )
        	sprintf( Pars.scene_psf_file, "%s%02d_psf.fits", Pars.root_name, position );
	else
		sprintf( Pars.scene_psf_file, "%s%03d_psf.fits", Pars.root_name, position );

	printf( "Reading input PSF from %s.\n", Pars.scene_psf_file );

        psf = Read_FITS( Pars.scene_psf_file, &nx, &ny );

        if ( nx != Pars.int_dim || ny != Pars.int_dim )
        {
                printf("ERROR : PSF dimensions do not match parameter file specs\n");
                exit(0);
        }

	Pars.current_pos = position;

	printf( "  Input critically-sampled undistorted PSF dimensions are %d by %d (%f arcsec/pixel).\n", nx, ny, Pars.psf_scale );

        return( psf );
}

/*--------------------------------------------------------------------------
* Read_scene :
*	Read user-generated undistorted ACS scene from a FITS file.
* Inputs/outputs :
*	nx, ny : pointers to integers that will return the image dimensions.
* This routine returns the image.  This routine allocates memory for the
* image using Alloc_image().
*--------------------------------------------------------------------------*/
float **Read_scene( int *nx, int *ny )
{
	float	**image, x_arcsec, y_arcsec;

	if ( strstr( Pars.scene_input_file, "blank_" ) == NULL )
	{
        	printf( "Reading input scene from %s.\n", Pars.scene_input_file );
        	image = Read_FITS( Pars.scene_input_file, nx, ny );
        	printf("  Input (undistorted) scene dimensions are %d by %d (%g arcsec/pixel).\n", 
			*nx, *ny, Pars.scene_pixel_scale );
	}
	else
	{
		printf( "Creating a blank input scene.\n" );
		sscanf( &Pars.scene_input_file[6], "%f_%f", &x_arcsec, &y_arcsec );
		*nx = ceil(x_arcsec / Pars.psf_scale);
		*ny = ceil(y_arcsec / Pars.psf_scale);
        	printf("  Input (undistorted) scene dimensions are %g by %g arcseconds.\n", x_arcsec, y_arcsec );
		image = NULL;
	}


	return( image );
}

/*----------------------------------------------------------------------------------
* Convolve_scene :
*	Resample the input scene, convolve it with the PSF, and map it to a
*	distorted grid.
* Inputs :
*	scene : Input scene
*	scene_nx, scene_ny : When called, contains the dimensions of the input scene
*	psf_in : PSF produced by tiny2
* Outputs :
*	scene_nx, scene_ny : Dimensions of processed image (previous values will
*			     be overwritten)
*	This routine returns the processed input file; memory is allocated with
*	Alloc_image();
*----------------------------------------------------------------------------------*/ 
float **Convolve_scene( float **scene, int *scene_nx, int *scene_ny, float **psf_in )
{
        int     num_grids = 4;
        int     grid_sizes[] = { 1024, 1280, 2048, 4096 };
        int     grid_dim, i, x, y, xx, yy, max_dim, new_scene_nx, new_scene_ny;
        float   c, mag, **convolved_scene, s_r, s_i, p_r, p_i;
        complex **scene_grid, **psf_grid;

	/* if input image is blank, don't resample or convolve */

	if ( strstr(Pars.scene_input_file, "blank_") != 0 )
	{
		convolved_scene = Alloc_image( *scene_nx, *scene_ny );
		return( convolved_scene );
	}

        /* determine if input scene needs to be resampled to match PSF sampling */

        if ( fabs((Pars.psf_scale-Pars.scene_pixel_scale)/Pars.psf_scale) > 0.005 )
        {
                mag = Pars.scene_pixel_scale / Pars.psf_scale;
                new_scene_nx = (int)(*scene_nx * mag + 1);
                new_scene_ny = (int)(*scene_ny * mag + 1);
        }
        else
        {
                mag = 1.0;
                new_scene_nx = *scene_nx;
                new_scene_ny = *scene_ny;
        }

        max_dim = new_scene_nx;
        if ( new_scene_ny > max_dim )
                max_dim = new_scene_ny;
        if ( Pars.int_dim > max_dim )
                max_dim = Pars.int_dim;

        /* Select convolution grid dimensions that can fully contain the PSF *
         * and the input scene; the dimensions are optimized for the FFT.    */

	grid_dim = grid_sizes[num_grids-1];

        for ( i = num_grids-1; i >= 0; --i )
                if ( grid_sizes[i] > max_dim )
                        grid_dim = grid_sizes[i];

	if ( new_scene_nx > grid_dim )
		new_scene_nx = grid_dim;
	if ( new_scene_ny > grid_dim )
		new_scene_ny = grid_dim;


        /* resample input scene, if necessary */

        if ( mag != 1.0 )
	{
		printf( "  Resampling input scene onto %d by %d grid (%f arcsec/pixel).\n", 
			grid_dim, grid_dim, Pars.psf_scale );
                scene_grid = Resample( scene, *scene_nx, *scene_ny, grid_dim, mag );
	}
        else
	{
		printf( "  Copying scene to %d by %d grid (%f arcsec/pixel).\n", 
			grid_dim, grid_dim, Pars.psf_scale );
                scene_grid = Copy_to_center( scene, *scene_nx, *scene_ny, grid_dim );
	}

        /* plop PSF into center of grid array */

	printf( "  Copying PSF to %d by %d grid (%f arcsec/pixel).\n", grid_dim, grid_dim, Pars.psf_scale );
        psf_grid = Copy_to_center( psf_in, Pars.int_dim, Pars.int_dim, grid_dim );


	printf( "  Computing FFT of scene.\n" );
	Shift_complex_to_origin( scene_grid, grid_dim );
        fft2d( scene_grid[0], grid_dim, -1 );

	printf( "  Computing FFT of PSF.\n" );
	Shift_complex_to_origin( psf_grid, grid_dim );
        fft2d( psf_grid[0], grid_dim, -1 );

	printf( "  Multiplying scene and PSF FFTs.\n" );

        for ( y = 0; y < grid_dim; ++y )
        {
                for ( x = 0; x < grid_dim; ++x )
                {
			s_r = scene_grid[y][x].r;
			s_i = scene_grid[y][x].i;
			p_r = psf_grid[y][x].r;
			p_i = psf_grid[y][x].i;

                        scene_grid[y][x].r = s_r*p_r - s_i*p_i;
                        scene_grid[y][x].i = s_i*p_r + s_r*p_i;
                }
        }

        Free_complex_image( psf_grid );

	printf( "  Taking inverse FFT.\n" );
        fft2d( scene_grid[0], grid_dim, 1 );
        Shift_complex_to_center( scene_grid, grid_dim );

        convolved_scene = Alloc_image( new_scene_nx, new_scene_ny );

	c = grid_dim * grid_dim;

        for ( y = 0; y < new_scene_ny; ++y )
        {
                yy = grid_dim/2 - new_scene_ny/2 + y;
                for ( x = 0; x < new_scene_nx; ++x )
                {
                        xx = grid_dim/2 - new_scene_nx/2 + x;
                        convolved_scene[y][x] = scene_grid[yy][xx].r / c;
                }
        }

        Free_complex_image( scene_grid );

        *scene_nx = new_scene_nx;
        *scene_ny = new_scene_ny;

	return( convolved_scene );
}

/*---------------------------------------------------------------------------------
*  Resample_thread :
*	Used by Resample; resamples a portion of an image within a thread
*---------------------------------------------------------------------------------*/
static void *Resample_thread( void *tpars )
{
	int	x, y;
	float	x_in, y_in;
	struct sparstruct *spars;
	struct isparstruct *ipars;

	spars = tpars;
	ipars = spars->ipars;

        for ( y = spars->start_row; y <= spars->end_row; ++y )
        {
                y_in = (y - ipars->dim_out/2) / ipars->mag + ipars->ny/2;
                if ( y_in < 0.0 || y_in >= ipars->ny )
                        continue;
                for ( x = 0; x < ipars->dim_out; ++x )
                {
                        x_in = (x - ipars->dim_out/2) / ipars->mag + ipars->nx/2;
                        if ( x_in < 0.0 || x_in >= ipars->nx )
                                continue;
                        ipars->newimage[y][x].r = Damped_interp( ipars->image, x_in, y_in, 
						ipars->nx, ipars->ny); /* , spars->thread ); */ 
                        ipars->newimage[y][x].i = 0.0;
                }
        }

	return( NULL );
}

/*---------------------------------------------------------------------------------
* Resample :
*	Resample an input image to a new grid.
* Inputs :
*	image : Image to be resampled
*	nx, ny : Dimensions of input image
*	dim_out : Dimensions of new, resampled image (dim_out by dim_out)
*	mag : magnification of newly resampled image
* Outputs :
*	This routine returns a new image (allocated with Alloc_image()) that
*	contains the resampled image.
*----------------------------------------------------------------------------------*/
complex **Resample( float **image, int nx, int ny, int dim_out, float mag )
{
	complex	**newimage;
	struct sparstruct spars[MAX_THREADS];
	struct isparstruct ispars;
	int	x, y;
	float	oldsum, sum;

#ifdef TT_THREADED
        pthread_t thread[MAX_THREADS];
        int     dy, ithread;
        void    *retval;
#endif

	oldsum = 0.0;
	for ( y = 0; y < ny; ++y )
		for ( x = 0; x < nx; ++x )
			oldsum += image[y][x];

	newimage = Alloc_complex_image( dim_out, dim_out );

	ispars.image = image;
	ispars.nx = nx;
	ispars.ny = ny;
	ispars.dim_out = dim_out;
	ispars.mag = mag;
	ispars.newimage = newimage;

#ifdef TT_THREADED
	dy = dim_out / Num_threads;

	for ( ithread = 0; ithread < Num_threads; ++ithread )
	{
		spars[ithread].ipars = &ispars;
		spars[ithread].thread = ithread;
		spars[ithread].start_row = ithread * dy;
                if ( ithread < Num_threads-1 )
                        spars[ithread].end_row = spars[ithread].start_row + dy - 1;
                else
                        spars[ithread].end_row = dim_out - 1;


                if ( pthread_create(&thread[ithread], NULL, Resample_thread, &spars[ithread]) )
                {
                        fprintf(stderr, "Cannot make thread %d\n", ithread);
                        exit(1);
                }
        }

        /* Join (collapse) the two threads */

        for ( ithread = 0; ithread < Num_threads; ++ithread )
        {
                if ( pthread_join(thread[ithread], &retval) )
                {
                        fprintf(stderr, "Thread join failed\n");
                        exit(1);
                }
        }
#else
        spars[0].ipars = &ispars;
        spars[0].thread = 0;
        spars[0].start_row = 0;
        spars[0].end_row = dim_out - 1;

        Resample_thread( &spars[0] );
#endif

	sum = 0.0;
	for ( y = 0; y < dim_out; ++y )
		for ( x = 0; x < dim_out; ++x )
			sum += newimage[y][x].r;
	for ( y = 0; y < dim_out; ++y )
		for ( x = 0; x < dim_out; ++x )
			newimage[y][x].r *= oldsum / sum;

	return( newimage );
}

/*---------------------------------------------------------------------------------
* Copy_to_center :
*	Create a blank image and copy another image into its center.
*--------------------------------------------------------------------------------*/
complex **Copy_to_center( float **image, int nx, int ny, int dim_out )
{
	complex **newimage;
	int	x, y, xx, yy;

	newimage = Alloc_complex_image( dim_out, dim_out );

	for ( y = 0; y < ny; ++y )
	{
		yy = dim_out/2 - ny/2 + y;

		for ( x = 0; x < nx; ++x )
		{
			xx = dim_out/2 - nx/2 + x;
			newimage[yy][xx].r = image[y][x];
			newimage[yy][xx].i = 0.0;
		}
	}

	return( newimage ); 
}

/*--------------------------------------------------------------------------------
* Read_optional_parameters :
*	Read the optional scene parameter file.
*-------------------------------------------------------------------------------*/
void Read_optional_parameters( void )
{
	FILE	*file;
	char	*line, keyword[MAX_STRING], value[MAX_STRING];
	int	iline, ip, xmax, ymax;

	strcpy( Pars.scene_input_file, "no_file" );
	Pars.scene_pixel_scale = 0.0;
	Pars.scene_x_center = -1;
	Pars.scene_y_center = -1;
	strcpy( Pars.scene_output_file, "no_file" );

	file = fopen( Pars.tt3_param_file, "r" );
	if ( file == NULL )
	{
		printf( "  Optional parameter file %s not found or cannot open.\n",
			Pars.tt3_param_file );
		exit(0);
	}

	ip = 0;
	iline = 0;
	while ( !feof(file) )
	{
		line = Get_string(file);
		if ( line == NULL )
			break;

		Get_keyword_value( line, keyword, value );
		if ( !strstr(keyword, "PSF") )
			strcompress( value );
		Str_up_case( keyword );

		if ( strstr(keyword, "INPUT_FILE") )
			strcpy( Pars.scene_input_file, value );
		else if ( strstr(keyword, "PIXEL_SCALE") )
			Pars.scene_pixel_scale = atof( value );
		else if ( strstr(keyword, "X_CENTER") )
			Pars.scene_x_center = atoi( value );
		else if ( strstr(keyword, "Y_CENTER") )
			Pars.scene_y_center = atoi( value );
		else if ( strstr(keyword, "OUTPUT_FILENAME") )
			strcpy( Pars.scene_output_file, value );
		else if ( strstr(keyword, "PSF") )
		{
			if ( ip >= MAX_PSFS )
			{
				printf( "ERROR : Number of PSFs in scene must be < %d\n", MAX_PSFS );
				exit(0);
			}

			sscanf( value, "%f %f %f", &Pars.scene_psf_x[ip], 
				&Pars.scene_psf_y[ip], &Pars.scene_psf_flux[ip] );
			++ip;
		}
		else
		{
			printf( "Read_optional_parameters : Unknown keyword in parameter file\n" );
			fclose( file );
			exit(0);
		}
		++iline;
	}

	fclose( file );

	Pars.num_scene_psfs = ip;

	if ( Pars.scene_pixel_scale <= 0.0 )
	{
		printf( "ERROR - Pixel scale in optional parameter file is <= 0.0\n" );
		exit(0);
	}

	if ( Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT || Pars.chip == ACS_SBC )
	{
		xmax = 1023;
		ymax = 1023;
	}
        else if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 )
        {
                xmax = 4095;
                ymax = 2050;
        }
        else if ( Pars.chip == WFC3_IR )
        {
                xmax = 1013;
                ymax = 1013;
        }
	else
	{
		xmax = 4095;
		ymax = 2047;
	}

	if ( Pars.scene_x_center < 0 || Pars.scene_x_center > xmax )
	{
		printf( "ERROR - Scene X position is outside of allowed range (0-%d)\n", xmax );
		exit(0);
	}

	if ( Pars.scene_y_center < 0 || Pars.scene_y_center > ymax )
	{
		printf( "ERROR - Scene Y position is outside of allowed range (0-%d)\n", ymax );
		exit(0);
	}


} /* Read_optional_parameters */

/*--------------------------------------------------------------------------------*/
void Add_psfs( float **psf, float **image, int nx, int ny )
{
	int	i, x_image_start, x_image_end, y_image_start, y_image_end;
	int	x_image, y_image;
	float	flux, xc_image, yc_image, x_psf_start, y_psf_start, x_psf, y_psf;

	/* psf and image scales are the same */

	/* psf is assumed to be flux normalized */

	for ( i = 0; i < Pars.num_scene_psfs; ++i )
	{
		flux = Pars.scene_psf_flux[i];
		xc_image = Pars.scene_psf_x[i] / Pars.psf_scale + nx/2;
		yc_image = Pars.scene_psf_y[i] / Pars.psf_scale + ny/2;

		if ( xc_image < 0 || roundval(xc_image) >= nx )
			continue;
		if ( yc_image < 0 || roundval(yc_image) >= ny )
			continue;

		x_image_start = roundval(xc_image - Pars.int_dim/2);
		x_image_end = x_image_start + Pars.int_dim - 1;
		if ( x_image_start < 0 )
			x_image_start = 0;
		if ( x_image_end >= nx )
			x_image_end = nx - 1;
		x_psf_start = Pars.int_dim/2 - (xc_image - x_image_start);

		y_image_start = roundval(yc_image - Pars.int_dim/2);
		y_image_end = y_image_start + Pars.int_dim - 1;
		if ( y_image_start < 0 )
			y_image_start = 0;
		if ( y_image_end >= ny )
			y_image_end = ny - 1;
		y_psf_start = Pars.int_dim/2 - (yc_image - y_image_start);

		y_psf = y_psf_start;
		for ( y_image = y_image_start; y_image <= y_image_end; ++y_image )
		{
			x_psf = x_psf_start;
			for ( x_image = x_image_start; x_image <= x_image_end; ++x_image )
			{
				image[y_image][x_image] += (flux * 
						Damped_interp( psf, x_psf, y_psf, Pars.int_dim, Pars.int_dim ));
				x_psf += 1.0;
			}
			y_psf += 1.0;
		}
	}

} /* Add_psfs */
