/* File     :  spectrum.c
 *
 * Contents :
 *
 * Author   :  John Krist (STScI)
 * Date     :  October 1998
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tinytim.h"

#define MAX_POINTS  15000
#define SQR(x) ((x)*(x))

#define DEFAULT_SPECTRUM 0
#define USER_SPECTRUM    1

#define FLUX_FLAM      1
#define FLUX_FNU       2
#define FLUX_JY        3
#define FLUX_PHOTLAM   4

/*---------------------------------------------------------------------------*/
static void Interpolate_spectrum( float *lambda, float *photlam, int npoints )
{
	int	i, j1, j2;
	float	dlambda, flux, slope, w, yint;

	/* lambda values must be increasing; spectrum fluxes in photons/Angstrom  *
	 * or other proportional units */

	j1 = 0;
	j2 = 0;
	for ( i = 0; i < Pars.num_waves; ++i )
	{
		w = Pars.wavelength[i];
		while ( lambda[j2] < w && j2 < npoints-1 )
			++j2;
		if ( j2 == 0 )
			j2 = 1;
		else
			j1 = j2 - 1;
		slope = (photlam[j2]-photlam[j1]) / (lambda[j2]-lambda[j1]);
		yint = photlam[j2] - slope * lambda[j2];
		flux = slope * Pars.wavelength[i] + yint;
		if ( flux < 0.0 )
			flux = 0.0;
		if ( i == 0 )
		  	dlambda = fabs(Pars.wavelength[1] - 
				       Pars.wavelength[0]);
		else if ( i == Pars.num_waves-1 )
		  	dlambda = fabs(Pars.wavelength[i] - 
				       Pars.wavelength[i-1]);
		else
		  	dlambda = fabs((Pars.wavelength[i-1] - 
				        Pars.wavelength[i+1])/2.0 );
		Pars.weight[i] *= flux * dlambda;
	}

} /* Interpolate_spectrum */

/*---------------------------------------------------------------------------*/
static void Interpol( float *lambda, float *photlam, int npoints )
{
	int	i, j1, j2;
	float	slope, w, yint;

	/* lambda values must be increasing; spectrum fluxes in photons/s/A */

	j1 = 0;
	j2 = 0;
	for ( i = 0; i < Pars.num_waves; ++i )
	{
		w = Pars.wavelength[i];
		while ( lambda[j2] < w && j2 < npoints-1 )
			++j2;
		if ( j2 == 0 )
			j2 = 1;
		else
			j1 = j2 - 1;
		slope = (photlam[j2]-photlam[j1]) / (lambda[j2]-lambda[j1]);
		yint = photlam[j2] - slope * lambda[j2]; 
                printf("interpol: slope=%f, wavelength=%f, yint=%f",slope,Pars.wavelength[i],yint);
		Pars.weight[i] *= (slope * Pars.wavelength[i] + yint);
	}

} /* Interpol */

/*---------------------------------------------------------------------------*/
static void Blackbody( void )
{
	int   i;
	double dlambda, flux, h, c, k, w, t, factor, exponent;
        int   useapproximation    = 0; 
        double const maxexpfactor = 50.0;

	t = 0.0;
	while ( t < 100.0 || t > 90000.0 )
	{	
		printf("\nEnter temperature (Kelvin) : ");                
		scanf("%lf", &t);
		if ( t < 100.0 || t > 90000.0 )
			printf("** Temperature must be 100 - 90000 K **\n"); 
	}

	h = 6.625e-27;
	c = 3.0e10;
	k = 1.38e-16;

        factor = 1.0; 

	for ( i = Pars.num_waves-1; i >=0; --i )
	{
		w = (double) Pars.wavelength[i] * 1.0e-4;  /* convert um to cm */
                
                /* checking the size of the exponent before to avoid floating point overflows */
                exponent = (h*c)/(w*k*t);
                
                if ((i==(Pars.num_waves-1)) && (exponent > maxexpfactor)) {
                      /* the exponent is too large to be computed directly. We can neglect the -1 and
                      compute the log of the function. We do normalize all subsequent fluxes to 
                      that of the first (largest) wave. */
                      factor = exponent;
                      useapproximation = 1;
                }

                if (useapproximation==1) {   
                   if ((-exponent + factor)< - maxexpfactor) {
                      /* the weight of this wavelength bin will be much smaller than
                      the precision of the averaging sum is anyway. In order to avoid
                      overflows, this flux can safely be set to 0.0, especially as
                      Pars.weight is float anyway.*/
                      flux = 0.0;
                   } else {
                      flux = exp(-log(pow(w,4.0)) - exponent + factor);
                   }
                } else {
		   flux = 1.0 / (pow(w,4.0) * (exp(exponent)-1.0));
                }   

		/* flux is now in phot/cm^2/s/A with arbitrary normalization */

		if ( i == 0 )
		  	dlambda = fabs(Pars.wavelength[1] - 
				       Pars.wavelength[0]);
		else if ( i == Pars.num_waves-1 )
		  	dlambda = fabs(Pars.wavelength[i] - 
				       Pars.wavelength[i-1]);
		else
		  	dlambda = fabs((Pars.wavelength[i-1] - 
				        Pars.wavelength[i+1])/2.0 );          
                                              
		Pars.weight[i] = (float)((double)Pars.weight[i] * flux * dlambda);                
	}

	sprintf( Pars.spectrum_file, "Blackbody(%fK)", t );

} /* Blackbody */

/*---------------------------------------------------------------------------
*  Power_law_nu :
*       Multiply the system throughput curve with a power law curve defined
*       as Flux(nu) = nu^alpha, where alpha is specified by the user.
*       The resulting fluxes are assumed to be directly proportional to
*       ergs/Hz (e.g. Janskys).  The normalization of the curve is not important.
*
*       The array Pars.weight[] must contain the system throughput
*       curve prior to calling this routine.
*--------------------------------------------------------------------------*/
static void Power_law_nu( void )
{
        double  alpha, c, dlambda, lambda, photlam, fnu, nu;
        int     i;


        printf( "\n\nYou have selected to use an object spectrum of the form : \n" );
        printf( "                   F(nu) = nu^alpha \n");
        printf( "where F(nu) is directly proportional to ergs/Hz.\n\n" );
        printf( "Enter the spectral index (alpha) : " );
        scanf( "%lf", &alpha );

        c = 3.0e14;     /* speed of light in microns */

        for ( i = 0; i < Pars.num_waves; ++i )
        {
                lambda = Pars.wavelength[i];
                nu = c / lambda;  /* convert microns to Hertz */
                fnu = pow(nu, alpha);

                /* "fnu" is currently in units ergs/Hz.  Convert to units directly      *
                 * proportional to photons/micron.  The normalization is not important. */

                photlam = fnu / lambda;

                if ( i == 0 )
                        dlambda = fabs(Pars.wavelength[1] -
                                       Pars.wavelength[0]);
                else if ( i == Pars.num_waves-1 )
                        dlambda = fabs(Pars.wavelength[i] -
                                       Pars.wavelength[i-1]);
                else
                        dlambda = fabs((Pars.wavelength[i-1] -
                                        Pars.wavelength[i+1])/2.0 );
                Pars.weight[i] *= photlam * dlambda;
        }
        sprintf( Pars.spectrum_file, "f(nu)=nu^%g", alpha );

} /* Power_law_nu */

/*---------------------------------------------------------------------------
*  Power_law_lambda :
*       Multiply the system throughput curve with a power law curve defined
*       as Flux(lambda) = lambda^beta, where beta is specified by the user.
*       The resulting fluxes are assumed to be directly proportional to
*       ergs/micron.  The normalization of the curve is not important.
*
*       The array Pars.weight[] must contain the system throughput
*       curve prior to calling this routine.
*--------------------------------------------------------------------------*/
static void Power_law_lambda( void )
{
        float   beta, dlambda, photlam, flam;
        int     i;


        printf( "\n\nYou have selected to use an object spectrum of the form : \n" );
        printf( "               F(lambda) = lambda^beta \n" );
        printf( "where F(lambda) is directly proportional to ergs/micron.\n\n");
        printf( "Enter the spectral index (beta) : " );
        scanf( "%f", &beta );

        for ( i = 0; i < Pars.num_waves; ++i )
        {
                flam = pow(Pars.wavelength[i], beta);

                /* flam is in ergs/lambda; convert to units directly  *
                 * proportional to photons/micron.  The normalization *
                 * is not important.                                  */

                photlam = flam * Pars.wavelength[i];

                if ( i == 0 )
                        dlambda = fabs(Pars.wavelength[1] -
                                       Pars.wavelength[0]);
                else if ( i == Pars.num_waves-1 )
                        dlambda = fabs(Pars.wavelength[i] -
                                       Pars.wavelength[i-1]);
                else
                        dlambda = fabs((Pars.wavelength[i-1] -
                                        Pars.wavelength[i+1])/2.0 );
                Pars.weight[i] *= photlam * dlambda;
        }

        sprintf( Pars.spectrum_file, "f(lambda)=lambda^%g", beta );

} /* Power_law_lambda */

/*---------------------------------------------------------------------------*/
static void Spectrum_table( int type )
{
	FILE 	*file;
	char 	line[MAX_STRING];
	int  	i, flux_type;
	float 	lambda[MAX_POINTS], photlam[MAX_POINTS], flux;
	double	c, h, photon_energy;


        c = 3.0e14;     /* speed of light in microns/sec */
        h = 6.6266e-27; /* Planck's constant */

	if ( type == DEFAULT_SPECTRUM )
	{
		Default_dir( Pars.spectrum_file, line );
		strcat( line, ".dat" );
		file = fopen( line, "r" );
		if ( file == NULL )
		{
			fprintf(stderr, "ERROR : default spectrum not found\n");
			exit(2);
		}
	}
	else
	{
		file = NULL;

		while ( file == NULL )
		{
			file = fopen( Pars.spectrum_file, "r" );
			if ( file == NULL )
			{
				fprintf(stderr, "\nERROR : Could not open %s\n", Pars.spectrum_file );
				printf("\nEnter name of spectrum file : ");
				scanf("%s", Pars.spectrum_file);
			}
		}
	}

        /* skip any comment lines (beginning with #) at beginning of the file */

        do {
                if ( fgets(line, MAX_STRING-1, file) == NULL )
                {
                        printf("ERROR : Error reading %s\n", Pars.spectrum_file);
                        fclose(file);
                        exit(0);
                }
        } while ( line[0] == '#' );

        /* First non-comment line in file must be FLAM, JY, FNU, or PHOTLAM */

        if ( strstr(line,"FLAM") != NULL )              /* ergs/micron */
                flux_type = FLUX_FLAM;
        else if ( strstr(line,"JY") != NULL )           /* Janskys */
                flux_type = FLUX_JY;
        else if ( strstr(line,"PHOTLAM") != NULL )      /* photons/micron */
                flux_type = FLUX_PHOTLAM;
        else if ( strstr(line,"FNU") != NULL )          /* ergs/Hz */
                flux_type = FLUX_FNU;
        else
        {
                printf( "ERROR : First non-comment line in spectrum file must be one\n");
                printf( "of the following : FNU, FLAM, JY, PHOTLAM\n");
                fclose( file );
                exit(0);
        }

        i = 0;
        while ( (fgets(line, MAX_STRING-1, file) != NULL) && (i < MAX_POINTS) )
        {
                /* must be wavelength(microns) & flux pairs */

                sscanf(line, "%f %f", &lambda[i], &flux);

		lambda[i] /= 10000.0;	/* convert Angstroms to microns */

                photon_energy = h * c / lambda[i];

                /* convert flux to photons/micron */

                if ( flux_type == FLUX_FLAM )
                        photlam[i] = flux / photon_energy;
                else if ( flux_type == FLUX_JY )
                        photlam[i] = 3.0e-9 * flux / SQR(lambda[i]) / photon_energy;
                else if ( flux_type == FLUX_FNU )
                        photlam[i] = 3.0e14 * flux / SQR(lambda[i]) / photon_energy;
                else
                        photlam[i] = flux;

                ++i;
        }

	fclose( file );

	if ( i >= MAX_POINTS )
	{
	   printf("\n ** ERROR : No. of points in spectrum (%d) > limit (%d) ** \n",
			i, MAX_POINTS );
	   exit(0);
	}

	Interpolate_spectrum( lambda, photlam, i );

} /* Spectrum_table */

/*---------------------------------------------------------------------------*/
static void Select_spectrum( void )
{
	char line[MAX_STRING];
	int  choice, i;
	char filename[MAX_SPECTRA][40];
	FILE *file;

	Default_dir( "spectra.lis", line );
	
	file = fopen( line, "r" );
	if ( file == NULL )
	{
	     fprintf(stderr, "Select_spectrum : Could not open spectra.lis\n");
	     exit(2);
	}

        fgets(line, MAX_STRING-1, file);
        fgets(line, MAX_STRING-1, file);

        printf("  #    Type    U-V    B-V    V-R    V-I    V-J    V-K\n");
        printf("-----------------------------------------------------\n");

        i = 0;
        while ( fgets(line, MAX_STRING-1, file) != NULL && i < MAX_SPECTRA )
        {
                sscanf(line, "%s", &filename[i][0]);
                printf(" %2d  %s", i+1, &line[9]);
                i = i + 1;
        }

        fclose( file );

	choice = 0;
	while ( choice < 1 || choice > i )
	{ 		
		printf("Enter spectrum # : ");
		scanf("%d", &choice);
	}

	strncpy(Pars.spectrum_file, &filename[choice-1][0], 39);

	Spectrum_table( DEFAULT_SPECTRUM );

} /* Select_spectrum */

/*---------------------------------------------------------------------------*/
static void Extinction( void )
{
	FILE 	*file;
	int	i, n;
	char	line[MAX_STRING], filename[MAX_STRING];
	float	c, lambda[50], ext[50];


	if ( Pars.av != 0.0 )
	{
		c = Pars.av;
		sprintf(line, "*(Av=%.2f mag)", c);
	}
	else
	{
		c = 3.1 * Pars.ebmv;
		sprintf(line, "*(E(B-v)=%.2f mag)", c);
	}
	strcat(Pars.spectrum_file, line);

	/* read in extinction curve from file */

	Default_dir( "extinct.dat", filename );

	file = fopen( filename, "r" );
	if ( file == NULL )
	{
		fprintf(stderr, "EXTINCTION: Could not open extinct.dat\n");
		exit(2);
	}

	fgets( line, MAX_STRING-1, file );
	while ( line[0] == '#' )
		fgets( line, MAX_STRING-1, file );
	sscanf( line, "%d", &n );	/* no. of points */
	fgets( line, MAX_STRING-1, file );
	while ( line[0] == '#' )
		fgets( line, MAX_STRING-1, file );

	sscanf( line, "%f %f", &lambda[0], &ext[0] );
	/* extinction in file is magnitudes; convert to flux ratio */
	ext[0] = pow( 10.0, -0.4*c*ext[0] );

	for ( i = 1; (i < n) && (i < 50); ++i )
	{
		fgets( line, MAX_STRING-1, file );
		sscanf( line, "%f %f", &lambda[i], &ext[i] );
		ext[i] = pow( 10.0, -0.4*c*ext[i] );
	}

	fclose( file );

	Interpol( lambda, ext, n );

} /* Extinction */

/*---------------------------------------------------------------------------*/
void Get_spectrum( void )
{
	int choice, i;
	double tot;

	choice = 0;
	while ( choice < 1 || choice > 5 )
	{
		printf("\nChoose form of object spectrum :\n");
	        printf("    1) Select a spectrum from list\n");
	        printf("    2) Blackbody\n");
                printf("    3) Power law : F(nu) = nu^i \n" );
                printf("    4) Power law : F(lambda) = lambda^i \n" );
	        printf("    5) Read user-provided spectrum from ASCII table\n");
        	printf("Choice : ");
		scanf( "%d", &choice );
	}

	switch ( choice )
	{
		case 1 : Select_spectrum();
			 break;
		case 2 : Blackbody();
			 break;
                case 3 : Power_law_nu();
                         break;
                case 4 : Power_law_lambda();
                         break;
		case 5 : printf("\nNOTE : Spectrum MUST have < %d points and be",
				MAX_POINTS );
			 printf("\n       in Angstrom and flux pairs (see manual for flux types)."); 
			 printf("\nEnter name of spectrum file : ");
			 scanf("%s", Pars.spectrum_file);
			 Spectrum_table( USER_SPECTRUM );
			 break;
	}

	/* if E(B-V) or Av are not zero, apply extinction */

	if ( Pars.ebmv != 0.0 || Pars.av != 0.0 )
	{
		printf("\nApplying extinction ...");
		Extinction();
	}

	/* normalize weights to total throughput of 1.0 */

	tot = 0.0;
	for ( i = 0; i < Pars.num_waves; ++i )
		tot = tot + (double) Pars.weight[i];
	for ( i = 0; i < Pars.num_waves; ++i )
		Pars.weight[i] /= (float) tot;

} /* Get_spectrum */

