/* File : distort.c
*
* Contents :
*	Damped_interp : Damped-sinc interpolation routine
*	Distort_thread : Threaded routine used by Distort
*	Distort : Map input image onto a distorted grid
*
* Written by John Krist, November 2000
*
* Typo in polynomial expression fixed, Richard Hook, March 2008
*
*/

#include <stdio.h>
#include <math.h>
#ifndef MACOSX
#include <malloc.h>
#endif
#include <stdlib.h>
#include "tinytim.h"

#ifdef TT_THREADED
#include <pthread.h>
#endif


struct idparstruct {
	float 	**in;
	int	nx, ny;
	int	x_det, y_det;
	float	**out;
	int	out_nx, out_ny;
	int	camera;
};
 
struct dparstruct {
	int	thread;
	int	start_row, end_row;
	struct idparstruct *ipars;
};

#define  K   13
#define  DK   6
#define  XK  21
#define  DXK 10

static float	Xk_threaded[MAX_THREADS][K], Yk_threaded[MAX_THREADS][K];
static float	Sinc_table[XK][XK][K][K];

/*--------------------------------------------------------------------------------------*/
void Init_damped_interp( void )
{
	int	i, ix, iy, j;
	float	xphase, yphase, sx, sy, dx, dy, c;


	for ( iy = 0; iy < XK; ++iy )
	{
		yphase = (iy - XK/2.0) / XK + 0.5/XK;
		for ( ix = 0; ix < XK; ++ix )
		{
			xphase = (ix - XK/2.0) / XK + 0.5/XK;

			for ( j = 0; j < K; ++j )
			{
				dy = j - DK - yphase;

				if ( fabs(dy) > DK )
				{
					sy = 0.0;
				}
				else if ( dy != 0.0 )
				{
					c = M_PI * dy;
					sy = sin(c)/c * sin(c/DK)/(c/DK);
				}
				else
				{
					sy = 1.0;
				}

				for ( i = 0; i < K; ++i )
				{
					dx = i - DK - xphase;

					if ( fabs(dx) > DK )
					{
						sx = 0.0;
					}
					else if ( dx != 0.0 )
					{
						c = M_PI * dx;
						sx = sin(c)/c * sin(c/DK)/(c/DK);
					}
					else
					{
						sx = 1.0;
					}

					Sinc_table[iy][ix][j][i] = sx * sy;
				}
			} 
		}
	}

} /* Init_damped_interp */

/*--------------------------------------------------------------------------------------*/
float Damped_interp( float **image, float x, float y, int nx, int ny )
{
	float	val;
	int	dx, dy, xk, yk, x0, x1, y0, y1, xpix, ypix;
	int	xa, ya, xb, yb, xi, yi;

	xpix = roundval(x);  /* coordinates are integral at center of pixel */
	ypix = roundval(y);
	if ( xpix < 0 || xpix >= nx || ypix < 0 || ypix >= ny )
		return( 0.0 );

	x0 = xpix - DK;
	y0 = ypix - DK;
	x1 = xpix + DK;
	y1 = ypix + DK;

	xa = 0;
	xb = K;
	if ( x0 < 0 )
		xa = -x0;
	if ( x1 >= nx )
		xb = nx - x1 + K - 1;

	ya = 0;
	yb = K;
	if ( y0 < 0 )
		ya = -y0;
	if ( y1 >= ny )
		yb = ny - y1 + K - 1;

	dx = (x - xpix + 0.5) * XK;
	if ( dx < 0 )
		dx = 0;
	else if ( dx >= XK )
		dx = XK - 1; 

	dy = (y - ypix + 0.5) * XK;
	if ( dy < 0 )
		dy = 0;
	else if ( dy >= XK )
		dy = XK - 1;

	val = 0.0;

	for ( yk = ya; yk < yb; ++yk )
	{
		yi = ypix - DK + yk;
		for ( xk = xa; xk < xb; ++xk )
		{
			xi = xpix - DK + xk;
			val += image[yi][xi] * Sinc_table[dy][dx][yk][xk];
		}
	}

	return( val );
}

/*--------------------------------------------------------------------------------------*/
float old_Damped_interp( float **image, float x, float y, int nx, int ny, int thread )
{
	float	dx, dy, t, val, c, *xk, *yk;
	int	i, ix, iy, x0, x1, y0, y1, xpix, ypix;
	int	xa, ya, xb, yb;

	xk = Xk_threaded[thread];
	yk = Yk_threaded[thread];

	xpix = roundval(x);
	ypix = roundval(y);
	if ( xpix < 0 || xpix >= nx || ypix < 0 || ypix >= ny )
		return( 0.0 );

	x0 = xpix - DK;
	y0 = ypix - DK;
	x1 = xpix + DK;
	y1 = ypix + DK;

	xa = 0;
	xb = K;
	if ( x0 < 0 )
		xa = -x0;
	if ( x1 >= nx )
		xb = nx - x1 + K - 1;

	ya = 0;
	yb = K;
	if ( y0 < 0 )
		ya = -y0;
	if ( y1 >= ny )
		yb = ny - y1 + K - 1;

	dx = xpix - x;
	dy = ypix - y;

	for ( i = 0; i < K; ++i )
	{
		t = i - DK + dx;
		if ( fabs(t) <= DK )
			if ( t != 0.0 )
			{
				c = M_PI * t;
				xk[i] = sin(c)/c * sin(c/DK)/(c/DK);
			}
			else
				xk[i] = 1.0;
		else
			xk[i] = 0.0;

		t = i - DK + dy;
		if ( fabs(t) <= DK )
			if ( t != 0.0 )
			{
				c = M_PI * t;
				yk[i] = sin(c)/c * sin(c/DK)/(c/DK);
			}
			else
				yk[i] = 1.0;
		else
			yk[i] = 0.0;
	}

	val = 0.0;
	for ( iy = ya; iy < yb; ++iy )
		for ( ix = xa; ix < xb; ++ix )
			val += image[y0+iy][x0+ix] * xk[ix] * yk[iy];

	return( val );

} /* old_Damped_interp */

/*--------------------------------------------------------------------------------------*/
static void *Distort_thread( void *tpars ) 
{
	float	y1, y2, xc, yc, xc0, yc0, xin, yin;
	double  xout, yout, dout;
	float   x, y, *cxc, *cyc, sample_area, subpixel_area=0.0;
	int	xout_pix, yout_pix, old_xout_pix, old_yout_pix, subfactor;
	int	center_in_pix_x, center_in_pix_y, center_out_pix_x, center_out_pix_y;
	struct	dparstruct *vpars;


	vpars = tpars;

	sample_area = Pars.psf_scale * Pars.psf_scale;

	subfactor = (int)(Pars.pixel_size / Pars.psf_scale * 3.0);
	subfactor = subfactor + 1 - subfactor % 2;
	if ( subfactor < 3 )
		subfactor = 3;
	dout = 1.0 / subfactor;


	center_in_pix_x = vpars->ipars->nx / 2;
	center_in_pix_y = vpars->ipars->ny / 2;
	center_out_pix_x = vpars->ipars->out_nx / 2;
	center_out_pix_y = vpars->ipars->out_ny / 2;

	cxc = Pars.xy_to_xc;
	cyc = Pars.xy_to_yc;

	x = vpars->ipars->x_det - Pars.x_ref;
	y = vpars->ipars->y_det - Pars.y_ref;
	xc0 =   cxc[0]*y + cxc[1]*x + 
		cxc[2]*y*y + cxc[3]*x*y + cxc[4]*x*x +
		cxc[5]*y*y*y + cxc[6]*x*y*y + cxc[7]*x*x*y + cxc[8]*x*x*x +
		cxc[9]*y*y*y*y + cxc[10]*x*y*y*y + cxc[11]*x*x*y*y + cxc[12]*x*x*x*y + cxc[13]*x*x*x*x;
	yc0 =   cyc[0]*y + cyc[1]*x + 
		cyc[2]*y*y + cyc[3]*x*y + cyc[4]*x*x +
		cyc[5]*y*y*y + cyc[6]*x*y*y + cyc[7]*x*x*y + cyc[8]*x*x*x +
		cyc[9]*y*y*y*y + cyc[10]*x*y*y*y + cyc[11]*x*x*y*y + cyc[12]*x*x*x*y + cyc[13]*x*x*x*x;

       /* Bug in above line fixed, RNH, March 2008 */

	old_xout_pix = -1;
	old_yout_pix = -1;

	y1 = vpars->start_row - 0.5 + dout/2;
	y2 = vpars->end_row + 0.5 - dout/2 + 0.01;

	for ( yout = y1; yout <= y2; yout += dout )
	{
		/* x,y = offset in pixels in distorted image from reference position */

		yout_pix = roundval(yout); 
		y = (yout - center_out_pix_y) / Pars.sub_factor + vpars->ipars->y_det - Pars.y_ref;

		for ( xout = dout/2-0.5; xout < vpars->ipars->out_nx+0.5-dout/2; xout += dout )
		{
			xout_pix = roundval(xout);

			if ( xout_pix != old_xout_pix || yout_pix != old_yout_pix )
			{
				subpixel_area = (dout * dout) * 
					Acs_pixel_area(xout_pix+vpars->ipars->x_det, yout_pix+vpars->ipars->y_det, 
							vpars->ipars->camera);
				old_xout_pix = xout_pix;
				old_yout_pix = yout_pix;
			}

			x = (xout - center_out_pix_x) / Pars.sub_factor + vpars->ipars->x_det - Pars.x_ref;


			xc =    cxc[0]*y + cxc[1]*x + 
				cxc[2]*y*y + cxc[3]*x*y + cxc[4]*x*x +
				cxc[5]*y*y*y + cxc[6]*x*y*y + cxc[7]*x*x*y + cxc[8]*x*x*x +
				cxc[9]*y*y*y*y + cxc[10]*x*y*y*y + cxc[11]*x*x*y*y + cxc[12]*x*x*x*y + cxc[13]*x*x*x*x;
			yc =    cyc[0]*y + cyc[1]*x + 
				cyc[2]*y*y + cyc[3]*x*y + cyc[4]*x*x +
				cyc[5]*y*y*y + cyc[6]*x*y*y + cyc[7]*x*x*y + cyc[8]*x*x*x +
				cyc[9]*y*y*y*y + cyc[10]*x*y*y*y + cyc[11]*x*x*x*y + cyc[12]*x*x*x*y + cyc[13]*x*x*x*x;


			xin = (xc - xc0) / Pars.psf_scale + center_in_pix_x;
			yin = (yc - yc0) / Pars.psf_scale + center_in_pix_y;

                        /*
			printf("-----------------------------------\n");
			printf("center_in_pix_x = %d\n", center_in_pix_x);
			printf("Pars.psf_scale = %f\n", Pars.psf_scale);
			printf("xout, yout = %f %f\n", xout, yout);
			printf("x, y = %f %f    xc, yc = %f %f\n", x, y, xc, yc);
			printf("xin, yin = %f %f\n", xin, yin);
                        */

			if ( xin >= 0 && xin < vpars->ipars->nx && yin >= 0 && yin < vpars->ipars->ny )
			{
			    if ( xout_pix < 0 || xout_pix >= vpars->ipars->out_nx || yout_pix < 0 || yout_pix >= vpars->ipars->out_ny )
				printf("xout_pix = %d  yout_pix = %d   nx = %d  ny = %d\n", xout_pix, yout_pix, 
					vpars->ipars->out_nx, vpars->ipars->out_ny );
			    else
				vpars->ipars->out[yout_pix][xout_pix] += 
					(Damped_interp(vpars->ipars->in, xin, yin, 
						vpars->ipars->nx, vpars->ipars->ny) / sample_area * subpixel_area);
			}
		}
	}

	return( NULL );

} /* Distort_thread */

/*--------------------------------------------------------------------------------------*/
void Distort( float **in, int nx, int ny, int x_det, int y_det, float **out, int out_nx, int out_ny, int camera )
{
	struct  idparstruct idpars;
        struct  dparstruct dpars[MAX_THREADS];

#ifdef TT_THREADED
        pthread_t thread[MAX_THREADS];
        int     dy, ithread;
        void    *retval;
#endif

	idpars.in = in;
	idpars.nx = nx;
	idpars.ny = ny;
	idpars.x_det = x_det;
	idpars.y_det = y_det;
	idpars.out = out;
	idpars.out_nx = out_nx;
	idpars.out_ny = out_ny;
	idpars.camera = camera;


#ifdef TT_THREADED
	dy = out_ny / Num_threads;

        for ( ithread = 0; ithread < Num_threads; ++ithread )
        {
                dpars[ithread].ipars = &idpars;
		dpars[ithread].thread = ithread;
                dpars[ithread].start_row = ithread * dy;
		if ( ithread < Num_threads-1 )
                	dpars[ithread].end_row = dpars[ithread].start_row + dy - 1;
		else
			dpars[ithread].end_row = out_ny - 1;


                if ( pthread_create(&thread[ithread], NULL, Distort_thread, &dpars[ithread]) )
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
	dpars[0].ipars = &idpars;
	dpars[0].thread = 0;
	dpars[0].start_row = 0;
	dpars[0].end_row = out_ny - 1;

	Distort_thread( &dpars[0] );
#endif

} /* Distort */

