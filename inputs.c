/* File      : inputs.c (formerly inputs2.c)
 *
 * Contents  :
 *	Find_string :  Search for a string in another string.
 *	Str_up_case  :  Convert a string to upper case.
 *	Strtoken   :  Get the next token from a string.
 *	Str_get_floats : Read an variable number of floats from
 *			a string.
 *	Open_filter_file : Open the filter data file.
 *	Get_filter : Get the wavelengths and throughput for a given filter.
 *	Read_positions : Read a position file.
 *	Read_table : Read wavelength/weight table.
 *
 * Author   : John Krist (STScI)
 * Date     : January 1992
 *
 * Modified:
 *
 *     R. N. Hook, March 1996
 *       Added NICMOS.
 *     J. Krist - December 1997
 *	 Added ACS, STIS
 *     J. Krist - October 1998
 *	 Replaced GetFilterwaves & GetFilterweights with GetFilter
 *     R. N. Hook & F. Stoehr, March 2008
 *       Add WFC3
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "tinytim.h"

#ifndef MACOSX
#include <malloc.h>
#endif
#include <memory.h>

/*-------------------------------------------------------------------------
*  Find_string : 
*       Find first occurrence of string2 in string1.  Returns pointer to
*       the location in string1 where string2 starts.  If string2 is 
*       not found, a NULL pointer is returned.    Check is case insensitive
*-------------------------------------------------------------------------*/
char *Find_string( char *string1a, char *string2a )
{
	char    *p0, *p1, *string1, *string2;

	string1 = (char *)malloc(strlen(string1a)+1);
	strcpy( string1, string1a );
	Str_up_case( string1 );
	string2 = (char *)malloc(strlen(string2a)+1);
	strcpy( string2, string2a );
	Str_up_case( string2 );

	p0 = string1;
	p1 = strstr( string1, string2 );
	free( string1 );
	free( string2 );

	if ( p1 == NULL )
		return( NULL );
	else
		return( &string1a[p1-p0] );
 
} /* Find_string */

/*--------------------------------------------------------------------------
*  Str_up_case :
*       Convert lowercase letters in string to upper case.
*--------------------------------------------------------------------------*/
void Str_up_case( char *string )
{
	while ( *string != '\0' )
	{
		*string = (char)toupper( *string );
		++string;
	}
}

/*------------------------------------------------------------------------
*  Str_token :
*	Return the next token in a string.  The token is placed in token
*	and a pointer to the character following that token in the input
*	string is returned.
*------------------------------------------------------------------------*/
char *Str_token( char *string, char *token )
{
	while ( isspace(*string) )
		++string;;

	while ( isgraph(*string) )
	{
		*token = *string;
		++token;
		++string;
	}

	*token = '\0';

	return( string );

} /* Str_token */

/*------------------------------------------------------------------------
*  Str_get_floats :
*	Read a variable number of floating point numbers from a string
*	into an array.
*------------------------------------------------------------------------*/
void Str_get_floats( char *string, int num_floats, float *floats )
{
	int	i;
	char	token[MAX_STRING];

	for ( i = 0; i < num_floats; ++i )
	{
		string = Str_token( string, token );
		sscanf( token, "%f", &floats[i] );
	}

} /* Str_get_floats */

/*------------------------------------------------------------------------
*  Open_filter_file :
*	Open the appropriate filter database file for the given camera.
*
*  Inputs :
*	chip  :  Camera number
*
*  Returns :
*	Returns a pointer to the opened file.
*------------------------------------------------------------------------*/
static FILE *Open_filter_file( int chip )
{
	char    filename[MAX_STRING]; 
	FILE    *file;

	if ( chip <= 8 )
		Default_dir( "wfpc1.fil", filename );
	else if ( chip >= WFPC2_PC && chip <= WFPC2_WFC_4 )
		Default_dir( "wfpc2.fil", filename );
	else switch( chip )
	{
	        case FOC_F48 : Default_dir( "focf48.fil", filename ); break;
	       case FOC2_F48 : Default_dir( "costarfocf48.fil", filename ); break;
	        case FOC_F96 : Default_dir( "focf96.fil", filename ); break;
	       case FOC2_F96 : Default_dir( "costarfocf96.fil", filename ); break;
	       case NICMOS_1 : Default_dir( "nicmos1.fil", filename ); break;
	       case NICMOS_2 : Default_dir( "nicmos2.fil", filename ); break;
	       case NICMOS_3 : Default_dir( "nicmos3.fil", filename ); break;
	    case NEWNICMOS_1 : Default_dir( "nicmos1.fil", filename ); break;
	    case NEWNICMOS_2 : Default_dir( "nicmos2.fil", filename ); break;
            case NEWNICMOS_3 : Default_dir( "nicmos3.fil", filename ); break;
	       case ACS_WFC1 : Default_dir( "acswfc.fil", filename ); break;
	       case ACS_WFC2 : Default_dir( "acswfc.fil", filename ); break;
	        case ACS_HRC : Default_dir( "acshrc.fil", filename ); break;
        case ACS_HRC_OFFSPOT : Default_dir( "acshrc.fil", filename ); break;
	        case ACS_SBC : Default_dir( "acssbc.fil", filename ); break;
	       case STIS_CCD : Default_dir( "stisccd.fil", filename ); break;
	       case STIS_NUV : Default_dir( "stisnuv.fil", filename ); break;
	       case STIS_FUV : Default_dir( "stisfuv.fil", filename ); break;
	     case WFC3_UVIS1 : Default_dir( "wfc3uvis.fil", filename ); break;
	     case WFC3_UVIS2 : Default_dir( "wfc3uvis.fil", filename ); break;
	        case WFC3_IR : Default_dir( "wfc3ir.fil", filename ); break;
	
		     default : printf("Error: Unknown camera mode\n"); exit(2);
	}

	if ( (file = fopen(filename, "r")) == NULL )
	{
		fprintf(stderr, "Open_filter_file : Error opening %s\n", filename);
		exit(2);
	}

	return( file );

} /* Open_filter_file */

/*------------------------------------------------------------------------
*  Interpolate_filter :
*	Resample filter curve via interpolation
*
*  Inputs :
*	old_num_waves :  Number of wavelengths in old curve
*	new_num_waves :  Number of wavelengths in new curve
*
*-------------------------------------------------------------------------*/
static void Interpolate_filter( int old_num_waves, int new_num_waves )
{
	float   old_weight[MAX_WAVES], old_wavelength[MAX_WAVES], w, slope, yint;
	float	dw;
	int	i, j1, j2;

	for ( i = 0; i < old_num_waves; ++i )
	{
		old_wavelength[i] = Pars.wavelength[i];
		old_weight[i] = Pars.weight[i];
	}

	dw = (old_wavelength[old_num_waves-1]-old_wavelength[0]) / (new_num_waves-1);

        /* note : wavelength values must be increasing */

        j1 = 0;
        j2 = 0;
	w = old_wavelength[0];

        for ( i = 0; i < new_num_waves; ++i )
        {
                while ( old_wavelength[j2] < w && j2 < old_num_waves-1 )
                        ++j2;
                if ( j2 == 0 )
                        j2 = 1;
                else
                        j1 = j2 - 1;
                slope = (old_weight[j2] - old_weight[j1]) / 
			(old_wavelength[j2] - old_wavelength[j1]);
                yint = old_weight[j2] - slope * old_wavelength[j2];
		Pars.wavelength[i] = w;
                Pars.weight[i] = slope * w + yint;
		w = w + dw;
        }

} /* Interpolate_filter */

/*------------------------------------------------------------------------
*  Get_filter :
*	Get the wavelengths and throughput for a camera/filter from database.
*
*  Inputs :
*
*  Returns :
*	Returns the number of wavelengths.
*-------------------------------------------------------------------------*/
int Get_filter( void )
{
	FILE	*file;
	char	line[MAX_STRING], *found_string;
	int	i, num_waves, new_num_waves;


	file = Open_filter_file( Pars.chip );

	/* Search for line containing filter name */

	do { 
		if ( fgets( line, MAX_STRING - 1, file ) == NULL )
		{
			fclose( file );
			return( 0 );
		}
		found_string = Find_string(line, Pars.filter_name);
	} while ( found_string == NULL );

	if ( found_string == NULL )
	{
		fclose( file );
		return( 0 );
	}

	sscanf( line, "%*s %d", &num_waves );  /* Get number of wavelengths */

	for ( i = 0; i < num_waves; ++i )
		fscanf( file, "%f %f", &Pars.wavelength[i], &Pars.weight[i] );

	fclose( file );

	if ( Pars.wave_mag != 1.0 )
	{
		new_num_waves = (float)num_waves * Pars.wave_mag;
		if ( new_num_waves > MAX_WAVES )
		{
		   fprintf(stderr, 
		      "\n** ERROR : Number of wavelengths in filter curve > 200\n");
		   exit(2);
		}
		Interpolate_filter( num_waves, new_num_waves );
		Pars.num_waves = new_num_waves;
	}
	else
		Pars.num_waves = num_waves;

	return( Pars.num_waves );

} /* Get_filter */

/*--------------------------------------------------------------------------
*  Read_positions :
*	Read X and Y positions from a file.  Return the number of positions
*	read.
*--------------------------------------------------------------------------*/
int Read_positions( char *filename, int *x, int *y )
{
	FILE    *file;
	int     i;
	char    line[MAX_STRING];
        float   xx, yy;

	if ( (file = fopen(&filename[1], "r")) == NULL )
	{
		fprintf(stderr, "Error opening %s\n", filename);
		exit(2);
	}

	i = 0;
	while( fgets(line, MAX_STRING-1, file) != NULL && i < MAX_POSITIONS )
	{
		sscanf(line, "%f %f", &xx, &yy);
                x[i] = (int) xx;
                y[i] = (int) yy;
		++i;
	}

	return( i );

} /* Read_positions */

/*--------------------------------------------------------------------------
*  Read_table :
*
* This module is intended to read a set of wavelength/weight pairs from
* a text file as a generalised way of creating a polychromatic PSF.
*
* The number of pairs read is returned.
*
* Richard Hook, October 1996
*
*--------------------------------------------------------------------------*/
int Read_table( char *filename, float *wave, float *weight )
{
	FILE    *file;
	int     i;
	char    line[MAX_STRING];

	if ( (file = fopen(filename, "r")) == NULL )
	{
		fprintf(stderr, "Error opening %s\n", filename);
		exit(2);
	}

	i = 0;
	while( fgets(line, MAX_STRING-1, file) != NULL && i < MAX_WAVES )
	{
		sscanf(line, "%f %f", &wave[i], &weight[i]);
		++i;
	}

	return( i );

} /* Read_table */


