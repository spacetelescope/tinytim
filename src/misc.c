/* File   :  misc.c
 *
 * Contents :
 *      round : round to nearest whole number
 *      Convolve : Convolve PSF with charge diffusion kernel
 *	Flux_normalize : flux normalizes an array
 *	sinc     :  Computes sinc(x)
 *      Str_trim :  Null terminate a string after last character.
 *      Swap_test:  Determine byte ordering of machine.
 *      strcompress : remove preceeding and trailing whitespace from string
 *      Get_keyword_value : return keyword and associated value from string
 *      Swap_int :  Byte swap an array of four byte integers/reals.
 *      Min_max  :  Compute min and max values in an image.
 *      xstrncpy :  Copy specified number of characters from one
 *                      array to another.
 *      Write_image :  Write a data array.
 *      Read_image  :  Read a data array created with Writeimage.
 *
 * Author   :  John Krist (STScI)
 * Date     :  January 1992
 *
 * Modifications :
 *	J. Krist - August 1994
 *		Modified image i/o to accept complex images.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifndef MACOSX
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include "tinytim.h"

/*-------------------------------------------------------------------------
*  roundval :
*       Round to nearest whole number.
*
*  This replaces rint() since rint appears not to be thread-safe under
*  Solaris (eventhough it is documented to be).
*-------------------------------------------------------------------------*/
int roundval( float value )
{
	return( (int)(value + 0.5) );
}

/*-------------------------------------------------------------------------
*  Convolve :
*	Convolve image with 3 x 3 kernel.
*
*  Inputs :
*	image : image array to be convolved.
*	kernel: 3 x 3 image array containing convolution kernel.
*	nx, ny : dimensions of image
*
*  Replaces image with convolved image.
*
*----------------------------------------------------------------------------*/
void Convolve( float **image, float kernel[3][3], int nx, int ny )
{
	float	*buffer[3], *temp, t, total;
	int	kern_dim, i, j, x, y;


	kern_dim = 3;  /* kernel is kern_dim x kern_dim */

	/* Allocate a 3 line buffer.  Since the convolution is done in place, *
	 * store needed, unconvolved lines in buffer.			      */

	for ( i = 0; i < kern_dim; ++i )
		buffer[i] = (float *)malloc( nx * sizeof(float) );

	fflush( stdout );

	/* Copy first three lines of image into buffers */

	for ( y = 0; y < kern_dim; ++y )
		for ( x = 0; x < nx; ++x )
			buffer[y][x] = image[y][x];

	/* Convolution does not occur at the edges of the array - those  *
	 * values are multiplied by the central kernel value.		 */

	for ( y = 1; y < ny-1; ++y )
	{
		for ( x = 1; x < nx-1; ++x )
		{
			total = 0.0;
			for ( j = 0; j < kern_dim; ++j )
				for ( i = 0; i < kern_dim; ++i )
				     total += buffer[j][x-1+i] * kernel[j][i];
			image[y][x] = total;
		}

		/* Shift lines in buffer and add new line */

		if ( y != ny-2 )
		{
			temp = buffer[0];
			buffer[0] = buffer[1];
			buffer[1] = buffer[2];
			buffer[2] = temp;
			for ( x = 0; x < nx; ++x )
				buffer[2][x] = image[y+2][x];
		}
	}

	t = kernel[kern_dim/2][kern_dim/2];

	for ( y = 1; y < ny-1; ++y )
	{
		image[y][0] *= t;
		image[y][nx-1] *= t;
	}

	for ( x = 0; x < nx-1; ++x )
	{
		image[0][x] *= t;
		image[ny-1][x] *= t;
	}

	for ( i = 0; i < kern_dim; ++i )
		free( buffer[i] );

} /* Convolve */

/*-------------------------------------------------------------------------*/
void Flux_normalize( float *array, int nx, int ny )
{
        long i, length;
        float sum;

        sum = 0.0;
        length = (long)nx * ny;

        for ( i = 0; i < length; ++i )
            sum += array[i];

        if ( sum != 0.0 )
            for ( i = 0; i < length; ++i )
                array[i] = array[i] / sum;
}

/*------------------------------------------------------------------*/
float sinc( float x )
{
        if ( x != 0.0 )
                return( sin(x) / x );
        else
                return( 1.0 );
} /* sinc */

/*----------------------------------------------------------------------*/
void Str_trim( char *string )
{
	/* Null terminate a string after the last alphanumeric character */

	char    *temp;

	if ( (temp = strchr( string, ' ' )) == NULL )
		return;
	else
		*temp = '\0';
}

/*----------------------------------------------------------------------*/
int Swap_test( void )
{
	/* Test the byte ordering of this machine.  Return 0 if most    *
	 * significant byte is first, or 1 if it is last (like on VAX).  *
	 * If ordering is unexpected (some wierd machine is being used), *
	 * return -1.                        */

	int     i;
	char    *temp;

	i = 1025;
	temp = (char *)&i;

	if ( temp[2] == 4 )
		return( 0 );
	else if ( temp[1] == 4 )
		return( 1 );
	else
		return( -1 );

} /* Swap_test */

/*-----------------------------------------------------------------------*/
void strcompress( char *string )
{
	int	i, j;

	i = 0;

	while ( isspace(string[i]) )
		++i;

	j = 0;
	while ( isgraph(string[i]) )
	{
		string[j] = string[i];
		++i;
		++j;
	}

	string[j] = '\0';
}

/*-----------------------------------------------------------------------*/
void Get_keyword_value( char *entry, char *keyword, char *value )
{
	size_t  eq_index;

	eq_index = strcspn( entry, "=" );

	strncpy( keyword, entry, eq_index );
	keyword[eq_index] = '\0';
	strcompress( keyword );

	strcpy( value, &entry[eq_index+1] );
}

/*-----------------------------------------------------------------------*/
void Swap_int( char *data, int number )
{
	char    temp;
	int     i;

	number = number * 4;    /* number of ints to number of bytes */

	for ( i = 0; i < number; i = i + 4 )
	{
		temp = data[i];
		data[i] = data[i+3];
		data[i+3] = temp;

		temp = data[i+1];
		data[i+1] = data[i+2];
		data[i+2] = temp;
	}
}

/*-----------------------------------------------------------------------*/
void Min_max( float *data, int num_pix, float *min_val, float *max_val )
{
	/* Compute the minimum and maximum values in the vector data,   *
	 * which has num_pix elements.                                   */

	int     i;

	*min_val = 1.0e32;
	*max_val = -1.0e32;

	for ( i = 0; i < num_pix; ++i )
	{
		if ( *data < *min_val )
			*min_val = *data;
		else if ( *data > *max_val )
			*max_val = *data;
		++data;
	}

} /* Min_max */

/*-----------------------------------------------------------------------*/
void xstrncpy( char *dest, char *source, int number )
{
	/* Copy number characters from source into dest, or until a      *
	 * \0 character is encountered.  If the number of characters     *
	 * in source is less than number, fill the rest of dest with     *
	 * spaces, up to number.  This does not copy the NULL character  *
	 * from source into dest.                                        */

	int     i;

	for ( i = 0; i < number; ++i )
		if ( source[i] != '\0' )
			dest[i] = source[i];
		else
			break;

	for ( i = i; i < number; ++i )
		dest[i] = ' ';

} /* xstrncpy */

/*-----------------------------------------------------------------------
*  Write_image :
*       Write an image to disk.
*
*  Inputs :
*       filename : Name of file to write out to.
*       image    : Array which contains the data.
*       nx, ny   : Dimensions of image.
*       datayype : Either FLOATING or COMPLEX.
*-----------------------------------------------------------------------*/
int Write_image( char *filename, void *image, int nx, int ny, int data_type )
{
	int     size;
	FILE    *file;

	if ( data_type == FLOATING )
		size = sizeof(float);
	else
		size = sizeof(complex);

	if ( (file = fopen( filename, "wb" )) == NULL )
		return( ERROR );

	if ( fwrite( image, nx*ny*size, 1, file ) != 1 )
	{
		fprintf(stderr, "Write_image - Error writing to file.\n");
		fclose( file );
		return( ERROR );
	}

	fclose( file );

	return( OKAY );

} /* Write_image */

/*-----------------------------------------------------------------------
*  Read_image :
*       Read an image written by Write_image into an array previously allocated
*       by Alloc_image or Alloc_complex_image.
*  Inputs :
*       filename : Name of file to read in.
*       image    : Array which will contain the data.
*       nx, ny   : Dimensions of image.
*       data_type : Either FLOATING or COMPLEX.
*-----------------------------------------------------------------------*/
int Read_image( char *filename, void *image, int nx, int ny, int data_type )
{
	FILE    *file;
	int     size;

	if ( (file = fopen( filename, "rb" )) == NULL )
		return( ERROR );

	if ( data_type == FLOATING )
		size = sizeof(float);
	else
		size = sizeof(complex);

	if ( fread( image, nx*ny*size, 1, file ) != 1 )
	{
		fclose( file );
		return( ERROR );
	}

	fclose( file );

	return( OKAY );

} /* Read_image */
