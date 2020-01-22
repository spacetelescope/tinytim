/* File     :  tiny1.c
 *
 * Contents :
 *	main()  :  Accept input parameters for Tiny Tim and write
 *		     them out to a file for use by tiny2.
 *
 * Author   :  John Krist (STScI)
 * Date     :  January 1992
 *
 * Modifications :
 *
 *	J. Krist - March 1992
 *	   Got rid of limit of 6.0 arcseconds for integrated PSF diameter.
 *	   Changed parameters for ComputeGridDims.
 *
 *      J. Krist - April 1992
 *	   Added WFPC II option and table of dates.
 *
 *      J. Krist - Aug 1992
 *	   Changed WFPC II WFC numbering from 1-3 to 2-4.
 *
 *      J. Krist - June 1993
 *	   Added some extra parameter default values.
 *
 *      J. Krist - March 1994
 *	   Modified some of the convergence criteria.
 *
 *      R. N. Hook - February 1996
 *         Modifications for NICMOS.
 *
 *      R. N. Hook - August 1996
 *         Set Pars.LimitRadius=3.0 for NICMOS, as recommended by John
 *         Krist to reduce SINC interpolation artifacts.
 *
 *      R. N. Hook - October 1996
 *         Added "TABLE" option to allow a list of wavelength/weight
 *         pairs to be read from a text file.
 *
 *      R. N. Hook - November 1996
 *         Incorporated several small enhancements suggested by John
 *         Krist. V4.2.
 *
 *	J. Krist - June 1997
 *	   Added query for X,Y position if NICMOS; added flags for
 *	   specifying whether NICMOS aberrations are computed using focus
 *	   and field position information
 *
 *	J. Krist - Sept 1997
 *	   Minor bug fixes
 *
 *	J. Krist - Nov 1997
 *	   Added the Advanced Camera, STIS
 *
 *	J. Krist - Oct 1998
 *	   Added inputs for new object spectrum feature;
 *	   added command line options for jitter and spectral extinction
 *
 *      J. Krist - Sept 2000
 *	   Altered ACS, STIS inputs; added WFC3
 *
 *      R. Hook & F. Stoehr - March 2008
 *         Added WFC3 support
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tinytim.h"

/*--------------------------------------------------------------------------
*  Process_command_line :
* 	Extract parameter file name and optional command line parameters
*
*-------------------------------------------------------------------------*/
static void Process_command_line( int argc, char *argv[], char *paramfile )
{
    int  	i, njitpars;
    float 	val;
    FILE 	*file;

    if ( argc < 2 )
    {
         printf("tiny1 paramfile [major=x minor=y angle=z] [jitter=x] [ebmv=x] [Av=x] [wmag=x]\n");
         exit(2);
    }

    strcpy( paramfile, argv[1] );

    if ( strcmp(paramfile, "-v") == 0 )
    {
        printf("%3.1f\n", (float)VERSION_NUM );
        exit(0);
    }

    Pars.jitter_major = -1.0;  /* default is no jitter (-1) */
    Pars.jitter_minor = -1.0;
    Pars.jitter_angle =  0.0;
    Pars.ebmv = 0.0;     /* default is no extinction */
    Pars.av = 0.0;
    Pars.wave_mag = 1.0;
    Pars.filter_file[0] = '\0';

    njitpars = 0;

    for ( i = 2; i < argc; ++i )
    {
        if ( Find_string(argv[i], "JITTER") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            if ( val < 0.5 || val > 200.0 )
            {
              fprintf(stderr,"\n** Jitter must be 0.5-200.0 mas\n");
              exit(2);
            }
            Pars.jitter_major = val;
            Pars.jitter_minor = val;
            Pars.jitter_angle = 0.0;
        }
        else if ( Find_string(argv[i], "MAJOR") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            if ( val < 0.5 || val > 200.0 )
            {
              fprintf(stderr,"\n** Jitter must be 0.5-200.0 mas\n");
              exit(2);
            }
            Pars.jitter_major = val;
            ++njitpars;
        }
        else if ( Find_string(argv[i], "MINOR") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            if ( val < 0.5 || val > 200.0 )
            {
              fprintf(stderr,"\n** Jitter must be 0.5-200.0 mas\n");
              exit(2);
            }
            Pars.jitter_minor = val;
            ++njitpars;
        }
        else if ( Find_string(argv[i], "ANGLE") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            Pars.jitter_angle = val;
            ++njitpars;
        }
        else if ( Find_string(argv[i], "EBMV") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            Pars.ebmv = val;
        }
        else if ( Find_string(argv[i], "AV") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            Pars.av = val;
        }
        else if ( Find_string(argv[i], "WMAG") != NULL )
        {
            sscanf( Find_string(argv[i], "="), "%*c%f", &val );
            Pars.wave_mag = val;
        }
        else if ( Find_string(argv[i], "FILTER") != NULL )
        {
            sscanf( Find_string(argv[i], "=")+1, "%s", Pars.filter_file );

            /* check and see if the filter file actually exists */

            file = fopen( Pars.filter_file, "r" );
            if ( file == NULL )
            {
                fprintf(stderr, "ERROR : Cannot open filter file %s\n",
                    Pars.filter_file );
                exit(2);
            }
            fclose( file );
        }
        else
        {
            fprintf(stderr, "\nUnrecognized parameter %s ignored\n",
                argv[i] );
        }
    }

    if ( njitpars != 0 && njitpars != 3 )
    {
        fprintf(stderr, "ERROR : Major and minor RMS jitter and major axis angle must all be\n");
        fprintf(stderr, "        specified for elliptical jitter.");
        exit(0);
    }
}

/*--------------------------------------------------------------------------
*  main() :
*	Accept inputs for Tiny Tim, compute grid sizes, and place the 
*	parameters into a file for use by tiny2.  This is the main routine
*	for tiny1.
*
*  Command line :
*	tiny1 paramfile [pupil=pupiltable] [aber=abertable]
*
*  where :
*	paramfile = Name of output parameter file (required).
*	pupiltable = Name of pupil table (optional; if not specified, then
*			the file TINYTIM/pupil.tab will be used.
*	abertable  = Name of aberration table (optional; if not specified,
*			then the file TINYTIM/aber.tab will be used.
*---------------------------------------------------------------------------*/
int main( argc, argv )
  int     argc;
  char    *argv[];
{
    struct  DateStruct obs_date;
    int     i, camera, oversampling_factor, new_dim;
    int	is_wfpc1, is_wfpc2, is_acs, is_nicmos, is_wfc3;
    char    position_file[MAX_STRING], yes_no[5];
        char    paramfile[MAX_STRING];
        char    manualfile[MAX_STRING];
    float   wavelength;
    float despace_in_microns, z4_offset_hack;

    is_wfpc1 = 0;
    is_wfpc2 = 0;
    is_nicmos = 0;
    is_acs = 0;
        is_wfc3 = 0;

    Process_command_line( argc, argv, paramfile );
    Default_dir( "tinytim.pdf", manualfile );

    printf("\n");
    printf("      ------ Tiny Tim v%3.1f : The HST PSF Generator ------\n",
        (float)VERSION_NUM);
    printf("\n");
    printf("                Release Date : June 1, 2008\n");
    printf("                   Developed by John Krist\n");
    printf("             Additional support by Richard Hook & Felix Stoehr\n");
    printf("        >> Manual in %s <<\n", manualfile);
    printf("        ** http://www.stsci.edu/software/tinytim **\n");

    camera = 0;
    while ( camera < 1 || camera > 23 )
    {
      printf("\n--------- Aberrated Cameras ---------  -------- Corrected HST Optics -------\n");
      printf("1) Wide Field Camera (WFPC1)             5) WFPC2 - Wide Field Camera \n");
      printf("2) Planetary Camera  (WFPC1)             6) WFPC2 - Planetary Camera \n");
      printf("3) FOC f/48                              7) COSTAR-corrected FOC f/48 \n");
      printf("4) FOC f/96                              8) COSTAR-corrected FOC f/96 \n");
      printf("\n--------- Second Generation ---------  -------- Third Generation ---------\n");
      printf(" 9) NICMOS Camera 1 (pre-cryocooler)    15) ACS - Wide Field Channel \n");
      printf("10) NICMOS Camera 2 (pre-cryocooler)    16) ACS - High Resolution Channel \n");
      printf("11) NICMOS Camera 3 (pre-cryocooler)    17) ACS - HRC coronagraph off-spot PSF\n");
      printf("12) STIS CCD                            18) ACS - Solar Blind Channel \n" );
      printf("13) STIS NUV-MAMA                       19) NICMOS Camera 1 + cryocooler \n");
      printf("14) STIS FUV-MAMA                       20) NICMOS Camera 2 + cryocooler \n");
      printf("                                        21) NICMOS Camera 3 + cryocooler \n");
      printf("--------- Fourth Generation --------- \n");
      printf("22) WFC3 UVIS channel\n");
      printf("23) WFC3 IR channel\n");
      printf("\nChoice : ");
      scanf("%d", &camera);
    }

    Pars.chip = 0;

    switch ( camera )
    {
         case 1 : while ( Pars.chip < 1 || Pars.chip > 4 )
             {
                printf("\nEnter chip (1,2,3,4) : ");
                 scanf("%d", &Pars.chip );
             }
             sprintf( Pars.camera_name, "WFC_chip_%d", Pars.chip);
             is_wfpc1 = 1;
             break;

        case 2 : while ( Pars.chip < 5 || Pars.chip > 8 )
              {
                printf("\nEnter chip (5,6,7,8) : ");
                 scanf("%d", &Pars.chip );
             }
             sprintf( Pars.camera_name, "PC_chip_%d", Pars.chip);
             is_wfpc1 = 1;
             break;

        case 3 : Pars.chip = FOC_F48;
             sprintf( Pars.camera_name, "FOC_f48" );
             break;

        case 4 : Pars.chip = FOC_F96;
             sprintf( Pars.camera_name, "FOC_f96" );
             break;

        case 5 : while ( Pars.chip < 2 || Pars.chip > 4 )
              {
                printf("\nEnter chip (2,3,4) : ");
                 scanf("%d", &Pars.chip );
             }
             sprintf( Pars.camera_name, "WFPC_II_WFC%d", Pars.chip);
             Pars.chip = Pars.chip - 2 + WFPC2_WFC_2;
             is_wfpc2 = 1;
             break;

        case 6 : Pars.chip = WFPC2_PC;
             sprintf( Pars.camera_name, "WFPC_II_PC" );
             is_wfpc2 = 1;
             break;

        case 7 : Pars.chip = FOC2_F48;
             sprintf( Pars.camera_name, "COSTAR+FOC_f48" );
             break;

        case 8 : Pars.chip = FOC2_F96;
             sprintf( Pars.camera_name, "COSTAR+FOC_f96" );
             break;

                case 9 : Pars.chip = NICMOS_1;
                         sprintf( Pars.camera_name, "NICMOS_1"  );
             is_nicmos = 1;
                         break;

                case 10 : Pars.chip = NICMOS_2;
                         sprintf( Pars.camera_name, "NICMOS_2" );
             is_nicmos = 1;
                         break;

                case 11 : Pars.chip = NICMOS_3;
                         sprintf( Pars.camera_name, "NICMOS_3" );
             is_nicmos = 1;
                         break;

        case 12 : Pars.chip = STIS_CCD;
             sprintf( Pars.camera_name, "STIS_CCD" );
             break;
        case 13 : Pars.chip = STIS_NUV;
             sprintf( Pars.camera_name, "STIS_NUV" );
             break;
        case 14 : Pars.chip = STIS_FUV;
             sprintf( Pars.camera_name, "STIS_FUV" );
             break;

           case 15 : while ( Pars.chip < 1 || Pars.chip > 2 )
               {
                printf("\nEnter detector (1 or 2) : ");
                 scanf("%d", &Pars.chip );
             }
             if ( Pars.chip == 1 )
             {
                 sprintf( Pars.camera_name, "ACS_WFC1" );
                Pars.chip = ACS_WFC1;
             }
             else
             {
                 sprintf( Pars.camera_name, "ACS_WFC2" );
                Pars.chip = ACS_WFC2;
             }
             is_acs = 1;
             break;
           case 16 : Pars.chip = ACS_HRC;
             sprintf( Pars.camera_name, "ACS_HRC" );
             is_acs = 1;
             break;
           case 17 : Pars.chip = ACS_HRC_OFFSPOT;
             sprintf( Pars.camera_name, "ACS_HRC_OFFSPOT" );
             is_acs = 1;
             break;
        case 18 : Pars.chip = ACS_SBC;
             sprintf( Pars.camera_name, "ACS_SBC" );
             is_acs = 1;
             break;

                case 19 : Pars.chip = NEWNICMOS_1;
                         sprintf( Pars.camera_name, "NEWNICMOS_1"  );
             is_nicmos = 1;
                         break;

                case 20 : Pars.chip = NEWNICMOS_2;
                         sprintf( Pars.camera_name, "NEWNICMOS_2" );
             is_nicmos = 1;
                         break;

                case 21 : Pars.chip = NEWNICMOS_3;
                         sprintf( Pars.camera_name, "NEWNICMOS_3" );
             is_nicmos = 1;
                         break;

        case 22 : while ( Pars.chip < 1 || Pars.chip > 2 )
                         {
                                printf("\nEnter detector (1 or 2) : ");
                                scanf("%d", &Pars.chip );
                         }
                         if ( Pars.chip == 1 )
                         {
                                sprintf( Pars.camera_name, "WFC3_UVIS1" );
                                Pars.chip = WFC3_UVIS1;
                         }
                         else
                         {
                                sprintf( Pars.camera_name, "WFC3_UVIS2" );
                                Pars.chip = WFC3_UVIS2;
                         }
                         is_wfc3 = 1;
                         break;

        case 23 : Pars.chip = WFC3_IR;
             sprintf( Pars.camera_name, "WFC3_IR" );
                         is_wfc3 = 1;
             break;
    }

    /* If camera is WFPC, WFPC2, NICMOS, ACS, WFC3 or STIS CCD get the  *
     * detector position at which to simulate the PSF.  A file     *
     * containing x,y pairs can be named instead by prefixing the  *
     * filename with a @.  If FOC, it is set to be at the center   *
     * of the detector (256,256).       			       */


/*        printf("%d %d %d %d %d",is_wfpc1, is_wfpc2,is_nicmos, is_acs ,is_wfc3) ;*/
    if ( is_wfpc1 || is_wfpc2 || is_nicmos || is_acs || is_wfc3 )
    {
               printf("\nEnter position (x and y) on detector in INTEGER\n");

        if ( is_wfpc1 || is_wfpc2 )
            printf("pixels (range = 0-799) or the filename of a list\n");
        else if ( is_nicmos )
               printf("pixels (range = 0-255) or the filename of a list\n");
        else if ( Pars.chip == ACS_WFC1 || Pars.chip == ACS_WFC2 )
        {
            printf("pixels (X range = 0-4095, Y range = 0-2047)\n");
            printf("or the filename of a list ");
        }
        else if ( Pars.chip == ACS_HRC || Pars.chip == ACS_HRC_OFFSPOT || Pars.chip == ACS_SBC )
        {
            printf("pixels (range = 0 to 1023) or the filename of a list\n" );
        }
                else if ( Pars.chip == WFC3_UVIS1 || Pars.chip == WFC3_UVIS2 )
                {
                        printf("pixels (X range = 0-4095, Y range = 0-2050)\n");
                        printf("or the filename of a list ");
                }
                else if ( Pars.chip == WFC3_IR )
                {
                        printf("pixels (range = 0 to 1013) or the filename of a list\n" );
                }
        else
            printf("pixels (range = 0 to 1023) or the filename of a list\n" );

        printf("of positions preceded by a '@' (ie. @xy.lis).\n");

        printf("\nPosition : ");
        scanf( "%s", position_file );

        if ( position_file[0] == '@' )
            Pars.num_pos = Read_positions( position_file, Pars.x, Pars.y );
        else
        {
            sscanf( position_file, "%d", &Pars.x[0]);
            scanf( "%d", &Pars.y[0]);
            Pars.num_pos = 1;
        }
    }
    else
    {
        Pars.x[0] = 0;
        Pars.y[0] = 0;
        Pars.num_pos = 1;
    }

    Pars.use_map = 2;

    Pars.old_nicmos = 0;
    if ( Pars.chip == NEWNICMOS_1 )
    {
        printf("\nIs the observation after May 16, 2002? (Y/N) : ");
        scanf( "%s", yes_no );
        if ( strchr(yes_no,'N') != NULL || strchr(yes_no,'n') != NULL )
            Pars.old_nicmos  = 1;
    }
    else if ( Pars.chip == NEWNICMOS_2 )
    {
        printf("\nIs the observation after September 29, 2002? (Y/N) : ");
        scanf( "%s", yes_no );
        if ( strchr(yes_no,'N') != NULL || strchr(yes_no,'n') != NULL )
            Pars.old_nicmos  = 1;
    }

    /* Get the camera parameters and aberration values */

    Init_entry_list();

    obs_date.day = 0;
    obs_date.month = 0;
    obs_date.year = 0;

    Read_pupil_table( Pars.chip,  Pars.pupil_name );
    Get_aberrations( &obs_date );

    /* Get the filter name, of the form Fxxx(W,LP,M,N).  If NONE is      *
     * given, then a monochromatic PSF will be generated at a wavelength *
     * specified by the user.                                            */

    Pars.num_waves = 0;

    if ( Pars.filter_file[0] != '\0' )
    {
        /* Read user-specified filter throughput table */

        printf("\nReading system throughput from %s", Pars.filter_file);
        printf("\nNOTE: Table should include ONLY the filter*system response.");
        printf("\n      Do not include object spectrum.  Table entries must be");
        printf("\n      wavelength (nm) and throughput pairs. \n");

            Pars.num_waves = Read_table( Pars.filter_file, Pars.wavelength, Pars.weight );
        strcpy( Pars.filter_name, Pars.filter_file );
        } 
    else
    {
        while ( Pars.num_waves == 0 )
        {
           printf("\nSelect filter passband :\n");
           printf("    - Enter the name of the filter (eg. f555w)\n");
           printf("    - Enter MONO to specify a single wavelength\n");
           printf("Filter : ");

            scanf("%s", Pars.filter_name);

           if ( Find_string(Pars.filter_name, "MONO") != NULL ||
                     Find_string(Pars.filter_name, "NONE") != NULL )
           {
            /* Compute a monochromatic PSF */

            Pars.num_waves = 1;
            Pars.weight[0] = 1.0;
            strcpy( Pars.spectrum_file, "No spectrum" );
                strcpy( Pars.filter_name, "No filter" );
            printf("\nCamera wavelength range is %.0f - %.0f nm\n",
                   Pars.min_wavelength, Pars.max_wavelength);

            wavelength = -1.0;

            while ( wavelength < 0 )
            {
               printf("Enter wavelength in nm : ");
               scanf( "%f", &wavelength );
               if ( wavelength < Pars.min_wavelength || wavelength > Pars.max_wavelength )
               {
               printf("\n** Error : wavelength is outside of ");
               printf("detector range (did you use NM?) **\n");
               wavelength = -1.0;
               }
               else
               Pars.wavelength[0] = wavelength;
            }
           }
           else
           {
             Pars.num_waves = Get_filter();

            if ( Pars.num_waves == 0 )
               printf("\n **ERROR : That filter is not in the database ** \n");
            else
            {
               printf("\n%d wavelengths will be used to generate PSF\n",
               Pars.num_waves );
            }
           }
             }
    }

    /* Convert wavelengths from nm to microns; check wavelength range */

    for ( i = 0; i < Pars.num_waves; ++i )
    {
            if ( Pars.wavelength[i] < Pars.min_wavelength ||
             Pars.wavelength[i] > Pars.max_wavelength )
        {
            printf("\n** WARNING : wavelength = %f nm is outside of camera's range",
                Pars.wavelength[i] );
            printf("\n             of %.0f - %.0f nm",
                Pars.min_wavelength, Pars.max_wavelength );
        }
        Pars.wavelength[i] /= 1000.0;
    }

    /* Select and apply object spectrum (if filter name was specified) */

    if ( Pars.num_waves > 1 )
        Get_spectrum();
    else
        strcpy( Pars.spectrum_file, "No spectrum" );

    /* If an ACS or WFC3 camera, tiny2 will compute a PSF that has a       *
     * sampling of the Nyquist sampling of the shortest wavelength *
     * in the specified filter passband, divided by 1.3.  Tiny3    *
     * will resample this onto distorted ACS pixels.	       */
 
    if ( is_acs || is_wfc3 )
        Pars.pixel_size = Pars.wavelength[0] * 360.0 * 3600.0 /
                    (4 * M_PI * DIAM_HST_MICRONS) / 1.3;

    /* Compute the Nyquist and integrated PSF grid dimensions */

    Compute_grid_dims( &Pars.psf_size, Pars.chip, Pars.num_waves,
        Pars.wavelength, Pars.crit_dim, &Pars.int_dim );

    /* Get the sampling factor (if not ACS or WFC3) */

    oversampling_factor = 1;

    if ( !is_acs && !is_wfc3 )
    {
        printf("\nDo you wish to generate a subsampled PSF? (Y/N) : ");
        scanf("%s", yes_no );

        if ( strchr(yes_no,'Y') != NULL || strchr(yes_no,'y') != NULL )
        {
              if ( is_wfpc2 || Pars.chip == STIS_CCD )
            {
                     printf("\n*** WARNING *** : The charge diffusion kernel is NOT\n");
                printf("applied to subsampled WFPC2, ACS, or STIS CCD PSFs.  You\n");
                printf("must convolve the PSF with the kernel (given in the FITS\n");
                printf("header) AFTER REBINNING TO NORMAL SAMPLING.\n");
            }

              oversampling_factor = 0;
              while ( oversampling_factor < 2 || oversampling_factor > 10 )
              {
                 printf("\nOverSampling factor along an axis (min=2, max=10) : ");
                scanf("%d", &oversampling_factor );
              }
        }
    }

    /* if normal sampling, apply charge diffusion kernel */

    Pars.scatter_flag = 1;

    if ( (is_wfpc2 || Pars.chip == STIS_CCD) && oversampling_factor != 1 )
        Pars.scatter_flag = 0;

    Pars.sampling = 1.0 / oversampling_factor;
    Pars.psf_scale = Pars.pixel_size / oversampling_factor;
    Pars.int_dim *= oversampling_factor;

    /* by default, compute field dependent aberrations */

    if ( is_wfpc2 || is_acs || is_wfc3 )
        Pars.adjust_for_xy = 1;

    /* if NICMOS, set flags so that aberrations are computed based  *
     * on focus and field position */

    if ( is_nicmos )
    {
        Pars.adjust_for_focus = 1;
        Pars.adjust_for_xy = 1;
    }

    if ( Pars.int_dim > 5e10 && !is_acs && !is_wfc3 )
    {
        printf("\nThe integrated PSF dimensions are %d by %d, which\n",
            Pars.int_dim, Pars.int_dim );
        printf("corresponds to %.2f arcsecs.\n", Pars.psf_size);
        printf("Do you want to decrease the PSF diameter? (Y/N) : ");
        scanf("%s", yes_no);
        if ( strchr(yes_no,'Y') != NULL || strchr(yes_no,'y') != NULL )
        {
            printf("Enter new PSF diameter in PIXELS (even #): ");
            scanf("%d", &new_dim);
            new_dim = new_dim - new_dim % 2;  /* force to be even */
            Pars.psf_size *= (float)new_dim / (float)Pars.int_dim;
            Pars.int_dim = new_dim;
        }
    }

    /* Add support for focus adustments, i.e. secondary mirror despace as
       described in the User's Guide, where:
              Z4 = OLD_Z4 + despace_in_microns * 0.011 
       and negative despace indicates movement of the secondary toward
       the primary.
     */
    printf("Secondary mirror despace is scaled by 0.011 and added to Z4.\n");
    printf("Focus, secondary mirror despace? [microns]: ");
    scanf("%f", &despace_in_microns);
    Hack_z4_entry(despace_in_microns * 0.011);

    printf("\nRootname of PSF image files (no extension) : ");
    scanf("%s", Pars.root_name);

    /* If SmartSkip set in parameter file, then the PSF for a given  *
     * wavelength will not be computed if its weight is less than    *
         * Pars.WeightLimit * max_weight			         */

    Pars.smart_skip = 1;
    Pars.weight_limit = 0.01;

    Write_params( paramfile, &obs_date );
    Delete_entry_list();

    if ( Pars.chip == NICMOS_3 )
        printf("\n** NOTE : Aberrations are for NICMOS 3 Campaign\n");

    printf("\nTo execute the next stage, enter the command :\n");
    printf("        tiny2 %s\n", paramfile);

    return(0);

} /* main */

