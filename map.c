/* File     :  map.c
 * 
 * Contents :
 *      Compute_map : Create the mirror zonal error map.
 *
 * Author   : John Krist
 * Date     : January 1992
 *
 * Modifications :
 *
 *	JEK - March 1994
 *		Added Wfpcmirrormap routine.
 *	JEK - April 1994
 *		Fixed bug in Wfpcmirrormap routine (did not zero arrays)
 *      RNH - March 1996
 *              Added NICMOS
 *      JEK - November 1996
 * 		Removed old maps.
 *      JEK - June 1997
 *              Modified NICMOS
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tinytim.h"

/*------------------------------------------------------------------------------
*  Compute_map :
*       This routine uses the new mirror maps.
*
*       Read in the primary and secondary mirror zonal error maps, resize them
*       to a specified dimension via bilinear interpolation, rotate them, and
*       add them together.
*
*  Inputs :
*       map_dim  :  Size of the mirror map array (same as the critically spaced
*                    PSF array).  Only the central map_dim/2 by map_dim/2 region 
*                    contains data.
*       chip    :  Camera id code.
*	theta   :  Rotation of camera.
*	v2, v3  :  Field position in arcsec in pupil coordinates
*  Returns:
*       Returns the map of mirror zonal errors in waves at 547 nm.
*---------------------------------------------------------------------------*/ 
float **Compute_map( int map_dim, int chip, float theta, float v2, float v3 )
{
	float   diam_illum_mm, **pri_map, **sec_map, **map1, **map2, t, xoff, yoff;
	float   t_v2, t_v3, pri_center, sec_center, map_center, new_radius, old_radius;
	float	v2_pc, v3_pc;
	int     x, y, pri_dim, sec_dim;
	char    map_name[MAX_STRING];

	/* secondary mirror map is aligned by default with the *
	 * WFPC2 PC camera center, so everything is relative   *
	 * to that position.				       */ 

	v2_pc = -2.374;
	v3_pc = 30.268;

	map_center = map_dim / 2.0;
	pri_dim = 344;			/* Primary is 280 x 280 pixels */
	pri_center = pri_dim / 2.0;
	sec_dim = 344;
	sec_center = sec_dim / 2.0;

	pri_map = Alloc_image( pri_dim, pri_dim );
	sec_map = Alloc_image( sec_dim, sec_dim );

	/* The mirror maps, as read in, are oriented so that +V3=right and   *
	 * +V2=up (in pupil coordinates), assuming the first pixel is in the *
	 * lower left.  The maps already account for reflection by           *
	 * multiplying the surface error by two.  They are in waves @        *
	 * 547 nm.				    			     */

	Default_dir( "primary.new", map_name );
	if ( Read_image( map_name, pri_map[0], pri_dim, pri_dim, FLOATING ) == ERROR )
	{
		fprintf(stderr, "Error reading file primary.new\n");
		exit(2);
	}

	Default_dir( "second.new", map_name );
	if ( Read_image( map_name, sec_map[0], sec_dim, sec_dim, FLOATING ) == ERROR )
	{
		fprintf(stderr, "Error reading file second.new\n");
		exit(2);
	}

	/* The primary map was obtained from phase retrieval, so the   *
	 * old_radius value is determined from the pupil magnification  *
	 * used in the fitting at 0.502 microns.		       */

	old_radius = 128.0 * 7.5 / (28.3 * 0.502 / 2.0);
	new_radius = map_dim / 4.0;	/* Radius of new map in pixels */

	/* The secondary map has already been interpolated to match the *
	 * spacing of the primary in pupil units.  Only a portion of    *
	 * the secondary is illuminated from a given field angle.       */

	diam_illum_mm = 267.3;            /* Diameter of beam at secondary in mm */

	/* t_v2 and t_v3 are shifts in mm on the secondary of the chief ray *
	 * from the OTA axis.  FocF96V3 and the like are in arcseconds.   *
	 * The value 4906.071 is the distance of the secondary from the   *
	 * primary in mm.						  */

	if ( (chip >= NICMOS_1 && chip <= NICMOS_3) || (chip >= NEWNICMOS_1 && chip <= NEWNICMOS_3) )
		   theta = theta + Pars.ota_offset;
	
	t_v2 = 4906.071 * tan((v2-v2_pc) * M_PI/(180.*3600.));
	t_v3 = 4906.071 * tan((v3-v3_pc) * M_PI/(180.*3600.));

	/* Convert to shifts in pixels in new map spacing */

	t_v2 = new_radius * t_v2 / (diam_illum_mm/2.0);
	t_v3 = new_radius * t_v3 / (diam_illum_mm/2.0);

	/* Scale and rotate primary map onto new grid */

	map1 = Rotate( pri_map, pri_dim, pri_dim, pri_center, pri_center, 
			map_dim, map_dim, map_center, map_center,
			theta, new_radius / old_radius );

	/* Scale, rotate, and shift secondary map onto new grid */

	t = -theta * M_PI / 180.0;
	xoff = t_v3 * cos(t) - t_v2 * sin(t);
	yoff = t_v3 * sin(t) + t_v2 * cos(t);

	map2 = Rotate( sec_map, sec_dim, sec_dim, sec_center, sec_center, 
			map_dim, map_dim, map_center - xoff, map_center - yoff,
			theta, new_radius / old_radius );

	/* Add the primary and secondary maps together. The maps are in  *
	 * waves at 547 nm.						 */

	for ( y = 0; y < map_dim; ++y )
		for ( x = 0; x < map_dim; ++x )
			map1[y][x] += map2[y][x];

	Free_image( pri_map );
	Free_image( sec_map );
	Free_image( map2 );

	return( map1 );

}  /* Compute_map */

