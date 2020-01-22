/* File : fitsio.c
 *
 * Contents :
 *	Which_byte_order : determines byte order of computer
 *	Report_error : prints error messages
 *	Read_FITS_header : read FITS file header and return image dimensions
 *	Swap_bytes : byte-swap 4-byte words
 *	Read_FITS : Read simple FITS file
 *	Write_FITS_header : Write basic image info to FITS header
 *	Write_FITS : Write image or PSF to FITS file
 *
 * Written by John Krist, November 2000
 * Updated for WFC3, Richard Hook & Felix Stoehr, March 2008
 */

#include <stdio.h>
#include <string.h>

#ifndef MACOSX
#include <malloc.h>
#include <stdlib.h>
#else
#include <stdlib.h>
#endif

#include <time.h>
#include "tinytim.h"

enum types { UNDEFINED, BIG_END, LIT_END };

/*-------------------------------------------------------------------------*/
static int Which_byte_order( void )
{
        /* Determines byte order of system.  If  *
         * system does not use IEEE (eg. VAX),   *
         * then the byte order is undetermined.  */

        float a = 1.4335;
        unsigned char big_end[] = {238, 124, 183, 63};
        unsigned char lit_end[] = {63, 183, 124, 238}; 
        unsigned char *b;
        int type;
        int   i;

        b = (unsigned char *)&a;
        type = UNDEFINED;

        if ( b[0] == big_end[0] )
        {
                for ( i = 1; i <= 3; ++i )
		{
                        if ( b[i] == big_end[i] )
                                type = BIG_END;
			else
			{
				type = UNDEFINED;
				break;
			}
		}
        }
        else if ( b[0] == lit_end[0] )
        {
                for ( i = 1; i <= 3; ++i )
		{
                        if ( b[i] == lit_end[i] )
                                type = LIT_END;
			else
			{
				type = UNDEFINED;
				break;
			}
		}
        }

        return( type );
}

/*-------------------------------------------------------------------------*/
static void Report_error( char *message, FILE *file )
{
  fprintf( stderr, message );
  fprintf( stderr, "\n" );
  if ( file != NULL )
    fclose( file );
  exit(-1);
}

/*-------------------------------------------------------------------------*/
static void Read_FITS_header( FILE *file, int *nx, int *ny )
{
	char	record[2880], *p;
	int	ival, line, nfound;

	/* NOTE : This routine is intended only to read in the header  *
	 * from the FITS files used by the Tiny Tim software.  It is   *
	 * not intended as a general FITS file header reader.          */

	if ( fread( record, 2880, 1, file ) != 1 )
		Report_error( "Error reading file header", file );

	nfound = 0;
	line = 0;
	p = record;

	while ( strncmp(p,"END ",4) != 0 ) 
	{
		if ( strncmp(p,"NAXIS",5) == 0 )
		{
			sscanf( p, "%*10c%d", &ival );
			if ( p[5] == ' ' )
			{
				if ( ival != 2 )
					Report_error("NAXIS not equal to 2", file);
				++nfound;
			}
			else if ( p[5] == '1' )
			{
				*nx = ival;
				++nfound;
			}
			else if ( p[5] == '2' )
			{
				*ny = ival;
				++nfound;
			}
		}
		else if ( strncmp( p, "BITPIX", 6 ) == 0 )
		{
			sscanf( p, "%*10c%d", &ival );
			if ( ival != -32 )
				Report_error( "BITPIX is not -32", file );
			++nfound;
		}
		p += 80;
		++line;

		if ( line == 36 )
		{
			if ( fread( record, 2880, 1, file ) != 1 )
				Report_error( "Error reading file header", file );
			line = 0;
			p = record;
		}
	}

	if ( nfound != 4 )
		Report_error( "Did not find a required header keyword", file );
}

/*---------------------------------------------------------------------------*/
static void Swap_bytes( float *image, long numpix )
{
	char	*p, temp;
	long	i;

	for ( i = 0; i < numpix; ++i )
	{
		p = (char *)&image[i];

		temp = p[0];
		p[0] = p[3];
		p[3] = temp;
		temp = p[1];
		p[1] = p[2];
		p[2] = temp;
	}
} 

/*-------------------------------------------------------------------------------
*  Read_FITS :
*	Read a FITS file produced by or supplied with Tiny Tim.
*	*** NOTE : This routine is intended to read in only Tiny Tim FITS
*	files; it is not designed to be a general FITS reader and will
*	fail completely on unexpected FITS file types (ie. doubles, etc.).
*
*  Inputs :
*	filename : string containing name of file to open (include extension)
*
*  Outputs :
*	  nx, ny : Dimension of image in FITS file
*
*  Returns :
*	2D image (memory is allocated by this routine;  caller is responsible
*	for freeing it via Free_image().
*------------------------------------------------------------------------------*/
float **Read_FITS( char *filename, int *nx, int *ny )
{
	FILE 	*file;
	long	numpix;
	float	**image;
	int	order;

	/* NOTE : This routine is intended to read in only those FITS files  *
	 * provided with the Tiny Tim software.  It is not meant to be a     *
	 * general FITS reader.						     */

	order = Which_byte_order();
	if ( order == UNDEFINED )
	{
		printf("ERROR : Undefined machine byte order\n");
		exit(0);
	}

	file = fopen( filename, "rb" );
	if ( file == NULL )
		Report_error( "Could not open FITS file", file );

	Read_FITS_header( file, nx, ny );

	numpix = (long)*nx * (long)*ny;

	image = Alloc_image( *nx, *ny );

	if ( fread( image[0], sizeof(float), numpix, file ) != numpix )
		Report_error( "Could not read all requested data", file );

	if ( order == BIG_END )
		Swap_bytes( image[0], numpix );

	fclose( file );

	return( image );

} /* Read_FITS */  

/*---------------------------------------------------------------------------*/
static void Header_string( char *keyword, char *value, char *comment, char *entry )
{
	int	c;

	memset( entry, ' ', 80 );

	sprintf( entry, "%-8s= '%s'", keyword, value );

	if ( comment != NULL )
	{
		c = 11 + strlen(value) + 2;
		entry[c] = ' ';
		if ( c < 31 )
			c = 31;
		else
			c = c + 1; 
		sprintf( &entry[c], "/ %s", comment );
	}
}

/*---------------------------------------------------------------------------*/
static void Write_FITS_header( FILE *file, int nx, int ny, int type )
{
	char	record[2880];
	int	camera, i, j;
	time_t	current_time;
	char	stime[80];

	memset( record, ' ', 2880 );

	sprintf( &record[80*0], "SIMPLE  =                    T" );
	sprintf( &record[80*1], "BITPIX  = %20d", -32 );
	sprintf( &record[80*2], "NAXIS   = %20d", 2 );
	sprintf( &record[80*3], "NAXIS1  = %20d", nx );
	sprintf( &record[80*4], "NAXIS2  = %20d", ny );
	sprintf( &record[80*5], "EXTEND  =                    F" );
	time( &current_time );
	strcpy( stime, ctime(&current_time) );
	stime[strlen(stime)-1] = '\0';

	Header_string( "CREATED", stime, "Time and date file was created", &record[80*6] );
	i = 6;

	if ( type == IMAGE_FILE )
	{
		Header_string( "OPT_FILE", Pars.tt3_param_file, "Optional parameter file used", &record[80*(i+1)] );
		Header_string( "IN_FILE",  Pars.scene_input_file, "Undistorted input scene", &record[80*(i+2)] );
		Header_string( "PSF_FILE", Pars.scene_psf_file, "Undistorted input PSF", &record[80*(i+3)] ); 
		sprintf( &record[80*(i+4)], "X_IMAGE = %20d / Image X center on detector", Pars.scene_x_center );
		sprintf( &record[80*(i+5)], "Y_IMAGE = %20d / Image Y center on detector", Pars.scene_y_center );
		i = i + 5;
	}

	if ( type == PSF_FILE || type == DISTORTED_PSF_FILE || type == IMAGE_FILE )
	{
		camera = Pars.chip;
		Header_string( "INSTRUME", Pars.camera_name, "Simulated instrument", &record[80*(i+1)] );

		sprintf( &record[80*(i+2)], "FOCUS   = %20.4f / PSF RMS focus (waves @ 547 nm)", Pars.focus );
		sprintf( &record[80*(i+3)], "X_COMA  = %20.4f / PSF RMS X-coma (waves @ 547 nm)", Pars.xcoma );
		sprintf( &record[80*(i+4)], "Y_COMA  = %20.4f / PSF RMS Y-coma (waves @ 547 nm)", Pars.ycoma );
		sprintf( &record[80*(i+5)], "X_ASTIG = %20.4f / PSF RMS 0d astig (waves @ 547 nm)", 
			Pars.xastig );
		sprintf( &record[80*(i+6)], "Y_ASTIG = %20.4f / PSF RMS 45d astig (waves @ 547 nm)", 
			Pars.yastig );
		sprintf( &record[80*(i+7)], "X_CLOVER= %20.4f / PSF RMS X-clover (waves @ 547 nm)", 
			Pars.xclover );
		sprintf( &record[80*(i+8)], "Y_CLOVER= %20.4f / PSF RMS Y-clover (waves @ 547 nm)", 
			Pars.yclover );
		sprintf( &record[80*(i+9)], "SPHERICL= %20.4f / PSF RMS spherical (waves @ 547 nm)", 
			Pars.spherical );

		if ( type == PSF_FILE )
			sprintf( &record[80*(i+10)], "PIXSCALE= %20.4f / Pixel scale in arcseconds", 
				Pars.psf_scale );
		else if ( type == DISTORTED_PSF_FILE )
			sprintf( &record[80*(i+10)], "PIXSCALE= %20.4f / Approximate pixel scale in arcseconds", 
				Pars.psf_scale );
		else
			i = i - 1;	/* do not print scale if image file */

		if ( Pars.num_waves > 1 )
		{
			Header_string( "SPECTRUM", Pars.spectrum_file, "PSF spectrum type or file", 
				&record[80*(i+11)] );
			Header_string( "FILTER", Pars.filter_name, "PSF filter", &record[80*(i+12)] );
			i = i + 12;
		}
		else
		{
			sprintf( &record[80*(i+11)],"WAVELNTH= %20f / PSF wavelength in microns",
				Pars.wavelength[0] );
			i = i + 11;
		}

		if ( (camera >= 1 && camera <= 8) || (camera >= WFPC2_PC && camera <= WFPC2_WFC_4) ||
		     (camera >= NICMOS_1 && camera <= NICMOS_3) || (camera >= ACS_WFC1 && camera <= ACS_SBC) ||
		     (camera == STIS_CCD) || (camera >= NEWNICMOS_1 && camera <= NEWNICMOS_3) )
		{
		    sprintf( &record[80*(i+1)], "X_PSF   = %20d / X position of PSF center in detector pixels",
			Pars.x[Pars.current_pos] );
		    sprintf( &record[80*(i+2)], "Y_PSF   = %20d / Y position of PSF center in detector pixels",
			Pars.y[Pars.current_pos] );
		    i = i + 2;
		}

		if ( Pars.num_comments != 0 )
		{
			for ( j = 0; j < Pars.num_comments; ++j )
			{
			     ++i;
			     sprintf( &record[80*i], "COMMENT   %-70.70s", Pars.comment[j] );
			}
		}
	}

	sprintf( &record[80*(i+1)], "END" );

	/* go through record and replace \0 with spaces */

	for ( i = 0; i < 2880; ++i )
		if ( record[i] == '\0' )
			record[i] = ' ';
	
	if ( fwrite( record, 1, 2880, file ) != 2880 )
		Report_error( "Could not write new header", file );  

} /* Write_FITS_header */

/*----------------------------------------------------------------------------------
*  Write_FITS :
*	Write simple FITS file.
*
*  Inputs :
*	filename : string containing filename (including extension)
*	image : floating point array containing image to write out
*	nx, ny : dimensions of image
*	type : image type (PSF, PUPIL, etc.)
*	overwrite : if routine needs to swap bytes and this flag is
*		    non-zero, then the swaps will occur in place,
*		    reordering the image; otherwise, a temporary
*		    image is created for the swapping
*----------------------------------------------------------------------------------*/ 
void Write_FITS( char *filename, float *image, int nx, int ny, int type, int overwrite )
{
	FILE	*file;
	long	i, numpix, remain;
	char	record[2880];
	int	order;
	float	*newimage, *temp;


	order = Which_byte_order();
	if ( order == UNDEFINED )
	{
		printf( "ERROR : Undefined byte order\n" );
		exit(0);
	}

	numpix = (long)nx * (long)ny;

	file = fopen( filename, "wb" );

	Write_FITS_header( file, nx, ny, type );

	newimage = NULL;

	if ( order == BIG_END )
	{
		if ( overwrite == 0 )
		{
			newimage = (float *)malloc( numpix * sizeof(float) );
			if ( newimage == NULL )
				Report_error( "Could not allocate temp image", file );

			for ( i = 0; i < numpix; ++i )
				newimage[i] = image[i];

			Swap_bytes( newimage, numpix );
			temp = image;
			image = newimage;
		}
		else
		{
			Swap_bytes( image, numpix );
		}
	}

	if ( fwrite( image, sizeof(float), numpix, file ) != numpix )
		Report_error( "Could not write data file", file );

	/* fill any unused bytes of the final 2880 byte record */

	/* remain = 2880 - (numpix * sizeof(float)) % 2880; */

	remain = 2880 - (numpix * sizeof(float)) % 2880;

	if ( remain != 2880 )
	{
		memset( record, 0, 2880 );
		if ( fwrite( record, 1, remain, file ) != remain )
			Report_error( "Could not finish writing data", file );
	}

	fclose( file );

	if ( newimage != NULL )
	{
		free( newimage );
		image = temp;
	}

} /* Write_FITS */
 
