/* File  :  makemaps.c
 *
 * Contents :
 *  	main()  :  Read the text files containing the primary and secondary 
 *		   mirror maps and write them out as binary files for use by 
 *		   Tiny Tim.
 *
 * Author : John Krist (STScI)
 * Date   : January 1992
 *
 * Modifications :
 *   J.E.K. - March 1994
 *	Added support for new wfpc mirror map.
 *
 *   J.E.K. - November 1996
 *      Removed old mirror maps.
 */

#include <stdio.h>
#include <stdlib.h>
#include "tinytim.h"

int main()
{
	FILE	*file;
	float	*map, **map0, constant;
	int	i, x, y;


	/* The mirror map data is in waves at 632.8 nm (historical reasons). *
         * It needs to be converted to waves at 547.0 nm.                    */

	constant = 632.8 / 547.0;

	/********************* primary.new ***********************/

	/* New primary map is 344 x 344 */

	map0 = Alloc_image( 344, 344 );
	map = map0[0];

	if ( (file = fopen("newpri.tab", "r")) == NULL )
	{
		fprintf(stderr, "Error opening file newpri.tab\n");
		exit(2);
	}

	printf("  Creating primary.new\n");

	for ( i = 0; i < 344*344; ++i )
		fscanf( file, "%f", map++ );

	fclose( file );

	for ( y = 0; y < 344; ++y )
		for ( x = 0; x < 344; ++x )
			map0[y][x] *= constant;

	i = Write_image( "primary.new", map0[0], 344, 344, FLOATING );

	Free_image( map0 );

	/********************* second.map ***********************/

	/* New secondary map is 344 x 344 */

	map0 = Alloc_image( 344, 344 );
	map = map0[0];

	if ( (file = fopen("newsec.tab", "r")) == NULL )
	{
		fprintf(stderr, "Error opening file newsec.tab\n");
		exit(2);
	}

	printf("  Creating second.new\n");

	for ( i = 0; i < 344*344; ++i )
		fscanf( file, "%f", map++ );

	fclose( file );

	for ( y = 0; y < 344; ++y )
		for ( x = 0; x < 344; ++x )
			map0[y][x] *= constant;

	i = Write_image( "second.new", map0[0], 344, 344, FLOATING );

	Free_image( map0 );

	return(0);
}
