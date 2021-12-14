#include "xstar.h"
#include "xstar_ext.h"


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com)
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/



/*
 * Gragg's method of Polynomial-Extrapolation with Modifications for
 * Conservation of Energy.  (i.e., this is a discrete mechanics method.)
 *
 * This method is from book _Numerical_Solutions_of_the_N-body_Problem_ by
 * Andrzej Marciniak.  This appears to be his favorite method.
 */

#define H	(1.)			/* step size			*/
#define P	(4)			/* half the order of the approx */
#define MAX_SKIP (32)			/* max # of correction skips	*/


/*
 * Notes:
 *
 * H should not be modified.  To change the step size, modify fv instead.
 *
 * P is one half the order of the approximation.  With P=4, FV=.25 you get
 * a very accurate approximation.  With P=2, FV=.25, you get an approximation
 * that is slightly worse than the taylor approximation with FV=2, but the
 * taylor method is much faster.  According to Andrzej's book, going beyond
 * P=4 won't get you much additional precision due to rounding errors.
 */


static int	ibeta[P+1];
static int	ibeta_2[P+1];
static double	beta[P+1];
static double	beta_2_inv[P+1];
static double	beta_inv[P+1];

static double	e0_2, a1, orig_a1, last_a1;
static double	alpha1[P+1], alpha2[P+1], bk[P+1], f[P], f1[P];
static point_2d	*x1, *y, *z, *txk[P+1];
static point_2d	*v1, *u, *w, *tvk[P+1];
static point_2d	*sum2, *sum3;
static double	*sum, *sum1;

static int	first_time = 1;

static int	eps_skip;
static int	max_eps_skip;


void restart_gpemce8( int n, double *m, point_2d *x, point_2d *v )
{
    return;
}

void init_gpemce8( int n, double *m, point_2d *x, point_2d *v )
{

    int		i, j, k, l;
    double	a, c;
    double	sum4, sum5;

    static int	last_n = -1;
    

    default_init( n, m, x, v );
    

    if( n != last_n )
    {
	last_n = n;
	
	/* allocate storage */
	if( !first_time )
	{
	    free( x1 );
	    free( y );
	    free( z );
	    free( v1 );
	    free( u );
	    free( w );
	    free( sum );
	    free( sum1 );
	    free( sum2 );
	    free( sum3 );
	    
	    for( k = 0; k < P+1; k++ )
	    {
		free( txk[k] );
		free( tvk[k] );
	    }
	}
	
	x1 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	y = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	z = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	v1 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	u = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	w = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	sum = (double *)Safe_Malloc( sizeof( double ) * n );
	sum1 = (double *)Safe_Malloc( sizeof( double ) * n );
	sum2 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	sum3 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	
	for( k = 0; k < P+1; k++ )
	{
	    txk[k] = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	    tvk[k] = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	}
    }
    

    /* calculate the total energy of the system (times 2) */
    e0_2 = 0;
    for( i = 0; i < n; i++ )
    {
	if( m[i] <= 0. )
	    continue;
	
	sum4 = 0;
	
	for( j = 0; j < n; j++ )
	{
	    if( m[j] <= 0. )
		continue;
	
	    if( i == j )
		continue;

	    
	    sum4 += m[j] / sqrt( (x[i].x - x[j].x)*(x[i].x - x[j].x)
				+ ( x[i].y - x[j].y )*( x[i].y - x[j].y ));
	}

	e0_2 += m[i]*( v[i].x*v[i].x + v[i].y*v[i].y - G*sum4 );
    }


    if( first_time )
    {
	/* step 1 */
	
	ibeta[0] = 1;
	ibeta[1] = 2;
	ibeta[2] = 3;
	for( k = 3; k < P+1; k++ )
	    ibeta[k] = 2*ibeta[k-2];

	for( k = 0; k < P+1; k++ )
	{
	    ibeta_2[k] = 2*ibeta[k];
	    beta[k] = ibeta[k];
	    beta_inv[k] = 1. / beta[k];
	    beta_2_inv[k] = .5 * beta_inv[k];
	}
	

	/* step 2 */
	for( k = 0; k < P+1; k++)
	{
	    bk[k] = 1;
	    
	    for( l = 0; l < P+1; l++ )
	    {
		if( k == l )
		    continue;
		
		bk[k] *= ibeta[l]*ibeta[l];
	    }
	}
	
	for( i = 0; i < P; i++ )
	    f[i] = f1[i] = 1;
	
	for( l = 0; l < P-1; l++ )
	{
	    for( i = 0; i < P-1; i++ )
	    {
		if( l == i )
		    continue;
		
		a = bk[l] - bk[i];
		f[l] *= (bk[P-1] - bk[i]) / a;
		f1[l] *= (bk[P] - bk[i]) / a;
	    }
	    
	    f[l] *= bk[P-1]/bk[l];
	    f1[l] *= bk[P]/bk[l];
	}
	
	a = c = 0;
	
	for( l = 0; l < P; l++ )
	{
	    a += f1[l];
	    c += f[l];
	}
	
	a = 1 - a;
	c = 1 - c;
	alpha1[P] = 0;
	alpha1[P-1] = 1/c;
	alpha2[P] = 1;
	alpha2[P-1] = -a/c;
	
	for( l = 0; l < P-1; l++ )
	{
	    alpha1[l] = -f[l]/c;
	    alpha2[l] = a*f[l]/c - f1[l];
	}
	
	sum4 = sum5 = 0;
	
	for( k = 0; k < P+1; k++ )
	{
	    a = pow( beta[k], 2*P );
	    sum4 += alpha1[k]/a;
	    sum5 += alpha2[k]/a;
	}
	
	orig_a1 = -sum4/sum5;
    }
    

    a1 = orig_a1;
    last_a1 = a1;
    eps_skip = -1;
    max_eps_skip = 4;
    
    first_time = 0;
	
}


void move_stars_gpemce8( point_2d *x, point_2d *x_new, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    double	eps = 1./(1024. * 1024. * 1024. * 1024. );
    
    int		ibetak_2;
    double	betak_2_inv;
    double	betak_inv;
    double	betak;

    int		p1;
    int		i, j, k;
    
    double	ap, b, c, df, ff, r, r1, r2;

    int		eps_calc;

    double	dx, dy;
    double	r_inv, r3_inv;

    double	m_i;

    int		xpt, ypt;
    
    point_2d	*v_new = vk[0];
    point_2d	*v = vk[1];
    
    

    /* step 3 */
    for( k = 0; k < P+1; k++ )
    {
	ibetak_2 = ibeta_2[k];
	betak = beta[k];
	betak_inv = beta_inv[k];
	betak_2_inv = beta_2_inv[k];
	
	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. )
		continue;
	    
	    /* step 4 */
	    z[i] = x[i];
	    u[i] = v[i];

	    /* step 5 */
	    y[i].x = z[i].x + H*betak_2_inv*u[i].x;
	    y[i].y = z[i].y + H*betak_2_inv*u[i].y;

	    sum2[i].x = sum2[i].y = 0;
	}

	/* step 5 (continued) */
	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. )
		continue;
	
	    for( j = i+1; j < max_stars; j++ )
	    {
		if( m[j] <= 0. )
		    continue;
		
		dx = z[i].x - z[j].x;
		dy = z[i].y - z[j].y;

		r3_inv = 1 / (dx*dx + dy*dy);
		r3_inv = sqrt( r3_inv ) * r3_inv;
		
		sum2[i].x += m[j] * dx * r3_inv;
		sum2[j].x -= m[i] * dx * r3_inv;
		sum2[i].y += m[j] * dy * r3_inv;
		sum2[j].y -= m[i] * dy * r3_inv;
	    }

	    w[i].x = u[i].x - H*G*betak_2_inv*sum2[i].x;
	    w[i].y = u[i].y - H*G*betak_2_inv*sum2[i].y;
	}
    
	
	for( p1 = 2; p1 < ibetak_2; p1++ )
	{
	    for( i = 0; i < max_stars; i++ )
	    {
		if( m[i] <= 0. )
		    continue;
		
		x1[i] = y[i];
		v1[i] = w[i];

		y[i].x = z[i].x + H*betak_inv*v1[i].x;
		y[i].y = z[i].y + H*betak_inv*v1[i].y;

		sum2[i].x = sum2[i].y = 0;
	    }
	    
	    for( i = 0; i < max_stars; i++ )
	    {
		if( m[i] <= 0. )
		    continue;
		
		for( j = i+1; j < max_stars; j++ )
		{
		    if( m[j] <= 0. )
			continue;
		
		    dx = x1[i].x - x1[j].x;
		    dy = x1[i].y - x1[j].y;
		    r3_inv = 1 / (dx*dx + dy*dy);
		    r3_inv = sqrt( r3_inv ) * r3_inv;
		
		    sum2[i].x += m[j] * dx * r3_inv;
		    sum2[j].x -= m[i] * dx * r3_inv;
		    sum2[i].y += m[j] * dy * r3_inv;
		    sum2[j].y -= m[i] * dy * r3_inv;
		}
		
		w[i].x = u[i].x - H*G*betak_inv*sum2[i].x;
		w[i].y = u[i].y - H*G*betak_inv*sum2[i].y;
	    }

	    for( i = 0; i < max_stars; i++ )
	    {
		z[i] = x1[i];
		u[i] = v1[i];
	    }
	}

	
	/* step 6 */
	for( i = 0; i < max_stars; i++ )
	{
	    sum2[i].x = sum2[i].y = 0;
	    sum3[i].x = sum3[i].y = 0;
	}

	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. )
		continue;
		
	    for( j = i+1; j < max_stars; j++ )
	    {
		if( m[j] <= 0. )
		    continue;
		
		dx = y[i].x - y[j].x;
		dy = y[i].y - y[j].y;
		r3_inv = 1 / (dx*dx + dy*dy);
		r3_inv = sqrt( r3_inv ) * r3_inv;
		
		sum2[i].x += m[j] * dx * r3_inv;
		sum2[j].x -= m[i] * dx * r3_inv;
		sum2[i].y += m[j] * dy * r3_inv;
		sum2[j].y -= m[i] * dy * r3_inv;
		
		dx = z[i].x - z[j].x + H*betak_inv*(w[i].x - w[j].x);
		dy = z[i].y - z[j].y + H*betak_inv*(w[i].y - w[j].y);
		r3_inv = 1 / (dx*dx + dy*dy);
		r3_inv = .5 * sqrt( r3_inv ) * r3_inv;
		
		sum3[i].x += m[j] * dx * r3_inv;
		sum3[j].x -= m[i] * dx * r3_inv;
		sum3[i].y += m[j] * dy * r3_inv;
		sum3[j].y -= m[i] * dy * r3_inv;
	    }

	    sum3[i].x += sum2[i].x;
	    sum3[i].y += sum2[i].y;

	    x1[i].x = (y[i].x + z[i].x)/2
		+ H*betak_2_inv*(w[i].x + u[i].x/2
				 - H*G*betak_2_inv*sum2[i].x );
	    x1[i].y = (y[i].y + z[i].y)/2
		+ H*betak_2_inv*(w[i].y + u[i].y/2
				 - H*G*betak_2_inv*sum2[i].y );

	    v1[i].x = (w[i].x + u[i].x)/2 - H*G*betak_2_inv*sum3[i].x;
	    v1[i].y = (w[i].y + u[i].y)/2 - H*G*betak_2_inv*sum3[i].y;
	}

	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. )
		continue;
	    
	    txk[k][i] = x1[i];
	    tvk[k][i] = v1[i];
	}
    }
    
    /* step 8 */
    for( eps_calc = 0; eps_skip <= 1 && eps_calc < 16; eps_calc++ )
    {
	ap = a1;

	ff = df = 0;
	
	for( i = 0; i < max_stars; i++ )
	    sum[i] = sum1[i] = 0;

	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. )
		continue;
	    
	    
	    for( j = i+1; j < max_stars; j++ )
	    {
		if( m[j] <= 0. )
		    continue;
		
		
		r1 = r2 = 0;
		for( k = 0; k < P+1; k++ )
		{
		    b = txk[k][i].x - txk[k][j].x;
		    r1 += (alpha1[k] + alpha2[k] * a1) * b;
		    r2 += alpha2[k] * b;
		}
		
		r = r1*r1;
		c = r1 * r2;
		
		r1 = r2 = 0;
		for( k = 0; k < P+1; k++ )
		{
		    b = txk[k][i].y - txk[k][j].y;
		    r1 += (alpha1[k] + alpha2[k] * a1) * b;
		    r2 += alpha2[k] * b;
		}
		
		r += r1*r1;
		c += r1 * r2;
		
		r3_inv = 1 / r;
		r_inv = sqrt( r3_inv );
		r3_inv *= c*r_inv;
		
		sum[i] += m[j] * r_inv;
		sum1[i] += m[j] * r3_inv;
		sum[j] += m[i] * r_inv;
		sum1[j] += m[i] * r3_inv;
	    }
	    
	    
	    r = c = 0;
	    
	    r1 = r2 = 0;
	    
	    for( k = 0; k < P+1; k++ )
	    {
		b = tvk[k][i].x;
		r1 += (alpha1[k] + alpha2[k] * a1) * b;
		r2 += alpha2[k] * b;
	    }
	    
	    r += r1*r1;
	    c += r1 * r2;
	    
	    r1 = r2 = 0;
	    
	    for( k = 0; k < P+1; k++ )
	    {
		b = tvk[k][i].y;
		r1 += (alpha1[k] + alpha2[k] * a1) * b;
		r2 += alpha2[k] * b;
	    }
	    
	    r += r1*r1;
	    c += r1 * r2;
	    
	    ff += m[i]*(r - G*sum[i]);
	    df += m[i]*(2*c + G*sum1[i]);
	}
	
	a1 -= (ff - e0_2) / df;

	if( fabs( a1 - ap ) < eps )
	    break;
    }
/*    if( eps_calc > 0 || fabs( a1 - last_a1 ) > .125*eps ) */
    if( eps_calc > 0 )
    {
	last_a1 = a1;
	if( eps_skip >= 0 )
	    max_eps_skip = (max_eps_skip + 1) >> 1;

	eps_skip = -2;
    }
    /*
     * we don't need to update a1 _every_ time, so skip a few.  updating a1
     * twice in a row is important.  It keeps the oscillations down.
     */
    eps_skip++;
    if( eps_skip > max_eps_skip )
    {
	eps_skip = 0;
	if( max_eps_skip < MAX_SKIP )
	    max_eps_skip += 1 + (max_eps_skip >> 3);
    }


    /* step 9 */
    for( i = 0; i < max_stars; i++ )
    {
	if( m[i] <= 0. || m[i] >= collapsar )
	{
	    x_new[i] = x[i];
	    v_new[i] = v[i];
	}
	else
	{
	    x_new[i].x = x_new[i].y = 0;
	    v_new[i].x = v_new[i].y = 0;
	}
    }

    for( k = 0; k < P+1; k++ )
    {
	b = alpha1[k] + alpha2[k] * a1;
	
	for( i = 0; i < max_stars; i++ )
	{
	    if( m[i] <= 0. || m[i] >= collapsar )
		continue;
		
	    x_new[i].x += b * txk[k][i].x;
	    x_new[i].y += b * txk[k][i].y;

	    v_new[i].x += b * tvk[k][i].x;
	    v_new[i].y += b * tvk[k][i].y;
	}
    }

    /*
     * plot the points
     */
    for( i = 0; i < max_stars; i++ )
    {
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;

	/* check for collisions */
	for( j = i+1; j < max_stars; j++ )
	{
	    if( m[j] <= 0.0 )
		continue;
	    
	    dx = x_new[i].x - x_new[j].x;
	    dy = x_new[i].y - x_new[j].y;
	    
	    if( dx*dx + dy*dy < collide_dist )
		collide( i, j, x, x_new, m, vk, ak, sp, sd );
	}
	
		
	xpt = x_new[i].x;
	ypt = x_new[i].y;
	
#	define cur	x_new
#	define prev	x
#	define v_0	v
#	include "plot_star.h"	
	
    }
}
