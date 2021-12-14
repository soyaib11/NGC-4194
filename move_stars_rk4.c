#include "xstar.h"
#include "xstar_ext.h"
#include <sys/time.h>


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
 * Runge-Kutta method of 4th order
 */

static point_2d	*ta_1, *ta_2, *ta_3;

static int	first_time = 1;

void restart_rk4( int n, double *m, point_2d *x, point_2d *v )
{
    return;
}


void init_rk4( int n, double *m, point_2d *x, point_2d *v )
{

    int		i;

    static int	last_n = -1;
    

    default_init( n, m, x, v );
    

    if( n != last_n )
    {
	last_n = n;
	
	/* allocate storage */
	if( !first_time )
	{
	    free( ta_1 );
	    free( ta_2 );
	    free( ta_3 );
	}
	
	ta_1 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	ta_2 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	ta_3 = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );

	for( i = 0; i < n; i++ )
	{
	    ta_1[i].x = ta_1[i].y = 0;
	    ta_2[i].x = ta_2[i].y = 0;
	    ta_3[i].x = ta_3[i].y = 0;
	}
    }

    first_time = 0;
	
}


#define collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd ) do_collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd, __FILE__, __LINE__ )

static void do_collide_rk4( int i, int j,
			   point_2d *cur, point_2d *prev,
			   double *m,
			   point_2d *vk[], point_2d *ak[],
			   point_2d *ta_1, point_2d *ta_2, point_2d *ta_3,
			   sys_param_type *sp, star_disp_type *sd,
			   char *file, int line
			   )
{
    /* collision:  totally inelastic -> combine */
    int		k;
    double tmass;
    double orig_m_i = m[i];
    
    num_collide++;
    if( verbose > 3 )
    {
	prt_sys_const( max_stars, prev, m, vk[1] );
	prt_sys_const_delta( orig_sc, max_stars, prev, m, vk[1] );
    }


    tmass = 1. / (m[i] + m[j]);
    
    for( k = 0; k < K; k++ )
    {
	vk[k][i].x = (m[i]*vk[k][i].x + m[j]*vk[k][j].x)*tmass;
	vk[k][i].y = (m[i]*vk[k][i].y + m[j]*vk[k][j].y)*tmass;
    }

    
    /* don't create a new collapsar */
    if( m[i] >= collapsar )
	m[i] += m[j];
    else if ( m[j] >= collapsar )
    {
	m[i] += m[j];

	cur[i] = cur[j];
	prev[i] = prev[j];
    }
    else
    {
	/* the new location of the combined stars */
	prev[i].x = (m[i]*prev[i].x + m[j]*prev[j].x)*tmass;
	prev[i].y = (m[i]*prev[i].y + m[j]*prev[j].y)*tmass;
	
	for( k = 0; k < K; k++ )
	{
	    ak[k][i].x = (m[i]*ak[k][i].x + m[j]*ak[k][j].x)*tmass;
	    ak[k][i].y = (m[i]*ak[k][i].y + m[j]*ak[k][j].y)*tmass;
	}
	ta_1[i].x = (m[i]*ta_1[i].x + m[j]*ta_1[j].x)*tmass;
	ta_1[i].y = (m[i]*ta_1[i].y + m[j]*ta_1[j].y)*tmass;
	ta_2[i].x = (m[i]*ta_2[i].x + m[j]*ta_2[j].x)*tmass;
	ta_2[i].y = (m[i]*ta_2[i].y + m[j]*ta_2[j].y)*tmass;
	ta_3[i].x = (m[i]*ta_3[i].x + m[j]*ta_3[j].x)*tmass;
	ta_3[i].y = (m[i]*ta_3[i].y + m[j]*ta_3[j].y)*tmass;
	
	m[i] += m[j];
	if( m[i] >= collapsar )
	{
	    /* This _should_, by construction, never happen. */
	    if( verbose > 0 )
		fprintf( stderr, "%s:%d  combined mass too large: m[%d]=%g\n",
			file, line, i, m[i]/fm );
	    
	    m[i] = collapsar * .999;
	}
    }
    
    --sp->live_stars;
    sd->tstamp = time(NULL);
    set_buffer_factor( sd, sp );
    if( verbose > 1 )
	fprintf( stderr, "%s:%d  live_stars=%d  stars: %d,%d  m: %g+%g=%g\n", file, line, sp->live_stars, i, j, orig_m_i/fm, m[j]/fm, m[i]/fm );
    
    m[j] = 0.0;
    cur[j].x = cur[j].y = prev[j].x = prev[j].y = -1.0;
    for( k = 0; k < K; k++ )
    {
	vk[k][j].x = vk[k][j].y = 0;
	ak[k][j].x = ak[k][j].y = 0;
    }	    
    ta_1[j].x = ta_1[j].y = 0;
    ta_2[j].x = ta_2[j].y = 0;
    ta_3[j].x = ta_3[j].y = 0;

    init_fna( max_stars, m, prev, vk[1] );
}



void move_stars_rk4( point_2d *prev, point_2d *cur, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    double		m_i, m_j;
    double		previ_x, previ_y;
    double		vi_x, vi_y;
    double		a0i_x, a0i_y;
    double		a1i_x, a1i_y;
    double		ai_x, ai_y;
    double		a_sum1, a_sum2;

    register int        i, j, k;	/* star index */

    int			xpt, ypt;

    point_2d	*v_0 = vk[0];
    point_2d	*v_1 = vk[1];
    point_2d	*a_0 = ak[0];
    
    

    /*
     * initialize the ta_[1-3] arrays
     */
    for( i = 0; i < max_stars; i++ )
    {
	ta_1[i].x = ta_1[i].y = 0;
	ta_2[i].x = ta_2[i].y = 0;
	ta_3[i].x = ta_3[i].y = 0;
    }
    


    MARK( main_loop );
    for( i = 0; i < max_stars; i++ )
    {
	MARK( main_loop2 );
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;
	
	previ_x = prev[i].x;
	previ_y = prev[i].y;

	/*
	 * calculate the a_0 values
	 */
	ai_x = a_0[i].x;
	ai_y = a_0[i].y;

	MARK( force_calc );
	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2, f;
	    register double f_d;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = prev[j].x - previ_x;
	    dy = prev[j].y - previ_y;
	    n2 = dx*dx + dy*dy;
	    
	    if( n2 < collide_dist )
	    {
		a_0[i].x = ai_x;
		a_0[i].y = ai_y;
		collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd );
		m_i = m[i];
		m_j = m[j];
		previ_x = prev[i].x;
		previ_y = prev[i].y;
		ai_x = a_0[i].x;
		ai_y = a_0[i].y;
	    }
	    else
	    {
		f = G / (sqrt(n2) * n2);
		
		f_d = dx * f;		
		ai_x += f_d * m_j;
		a_0[j].x -= f_d * m_i;
		
		f_d = dy * f;
		ai_y += f_d * m_j;
		a_0[j].y -= f_d * m_i;
	    }
	}
	a_0[i].x = ai_x;
	a_0[i].y = ai_y;
	
	
	/*
	 * calculate the ta_1 values
	 */
	ai_x = ta_1[i].x;
	ai_y = ta_1[i].y;
	vi_x = v_1[i].x;
	vi_y = v_1[i].y;

	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2, f;
	    register double f_d;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = prev[j].x - previ_x + .5*(v_1[j].x - vi_x);
	    dy = prev[j].y - previ_y + .5*(v_1[j].y - vi_y);

	    n2 = dx*dx + dy*dy;
	    
	    if( n2 < collide_dist )
	    {
		ta_1[i].x = ai_x;
		ta_1[i].y = ai_y;
		collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd );
		m_i = m[i];
		m_j = m[j];
		previ_x = prev[i].x;
		previ_y = prev[i].y;
		vi_x = v_1[i].x;
		vi_y = v_1[i].y;
		ai_x = ta_1[i].x;
		ai_y = ta_1[i].y;
	    }
	    else
	    {
		f = G / (sqrt(n2) * n2);
		
		f_d = dx * f;		
		ai_x += f_d * m_j;
		ta_1[j].x -= f_d * m_i;
		
		f_d = dy * f;
		ai_y += f_d * m_j;
		ta_1[j].y -= f_d * m_i;
	    }
	}
	ta_1[i].x = ai_x;
	ta_1[i].y = ai_y;
    }
    
    
    MARK( main_loop_b );
    for( i = 0; i < max_stars; i++ )
    {
	MARK( main_loop2_b );
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;
	
	previ_x = prev[i].x;
	previ_y = prev[i].y;

	/*
	 * calculate the ta_2 values
	 */
	vi_x = v_1[i].x;
	vi_y = v_1[i].y;
	a0i_x = a_0[i].x;
	a0i_y = a_0[i].y;
	ai_x = ta_2[i].x;
	ai_y = ta_2[i].y;
	
	MARK( force_calc_b );
	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2, f;
	    register double f_d;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = prev[j].x - previ_x + .5*(v_1[j].x - vi_x)
		+ .25*(a_0[j].x - a0i_x);
	    dy = prev[j].y - previ_y + .5*(v_1[j].y - vi_y)
		+ .25*(a_0[j].y - a0i_y);

	    n2 = dx*dx + dy*dy;
	    
	    if( n2 < collide_dist )
	    {
		ta_2[i].x = ai_x;
		ta_2[i].y = ai_y;
		collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd );
		m_i = m[i];
		m_j = m[j];
		previ_x = prev[i].x;
		previ_y = prev[i].y;
		vi_x = v_1[i].x;
		vi_y = v_1[i].y;
		a0i_x = a_0[i].x;
		a0i_y = a_0[i].y;
		ai_x = ta_2[i].x;
		ai_y = ta_2[i].y;
	    }
	    else
	    {
		f = G / (sqrt(n2) * n2);
		
		f_d = dx * f;		
		ai_x += f_d * m_j;
		ta_2[j].x -= f_d * m_i;
		
		f_d = dy * f;
		ai_y += f_d * m_j;
		ta_2[j].y -= f_d * m_i;
	    }
	}
	ta_2[i].x = ai_x;
	ta_2[i].y = ai_y;


	/*
	 * calculate the ta_3 values
	 */
	a1i_x = ta_1[i].x;
	a1i_y = ta_1[i].y;
	ai_x = ta_3[i].x;
	ai_y = ta_3[i].y;
	
	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2, f;
	    register double f_d;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = prev[j].x - previ_x + v_1[j].x - vi_x + .5*(ta_1[j].x - a1i_x);
	    dy = prev[j].y - previ_y + v_1[j].y - vi_y + .5*(ta_1[j].y - a1i_y);

	    n2 = dx*dx + dy*dy;
	    
	    if( n2 < collide_dist )
	    {
		ta_3[i].x = ai_x;
		ta_3[i].y = ai_y;
		collide_rk4( i, j, cur, prev, m, vk, ak, ta_1, ta_2, ta_3, sp, sd );
		m_i = m[i];
		m_j = m[j];
		previ_x = prev[i].x;
		previ_y = prev[i].y;
		vi_x = v_1[i].x;
		vi_y = v_1[i].y;
		a1i_x = ta_1[i].x;
		a1i_y = ta_1[i].y;
		ai_x = ta_3[i].x;
		ai_y = ta_3[i].y;
	    }
	    else
	    {
		f = G / (sqrt(n2) * n2);
		
		f_d = dx * f;		
		ai_x += f_d * m_j;
		ta_3[j].x -= f_d * m_i;
		
		f_d = dy * f;
		ai_y += f_d * m_j;
		ta_3[j].y -= f_d * m_i;
	    }
	}
	ta_3[i].x = ai_x;
	ta_3[i].y = ai_y;
	
	
	/*
	 * Move the star
	 */
	
	MARK( move_star );
	if( m_i >= collapsar )
	{
	    /* I think it looks better if the collapsars don't move ;-) */
	    v_1[i].x = v_1[i].y = 0.0;
	    
	    a_0[i].x = a_0[i].y = 0.0;
	    ai_x = ai_y = 0.0;
	    ta_1[i].x = ta_1[i].y = 0.0;
	    
	    xpt = cur[i].x = prev[i].x;
	    ypt = cur[i].y = prev[i].y;
	}
	else
	{
	    a_sum1 = ta_1[i].x + ta_2[i].x;
	    a_sum2 = (1./6)*(a_0[i].x + a_sum1);
	    xpt = cur[i].x = previ_x + (v_1[i].x + a_sum2);
	    v_0[i].x = v_1[i].x + a_sum2 + (1./6)*(a_sum1 + ai_x);
	    
	    a_sum1 = ta_1[i].y + ta_2[i].y;
	    a_sum2 = (1./6)*(a_0[i].y + a_sum1);
	    ypt = cur[i].y = previ_y + (v_1[i].y + a_sum2);
	    v_0[i].y = v_1[i].y + a_sum2 + (1./6)*(a_sum1 + ai_y);
	}
	ak[K-1][i].x = ak[K-1][i].y = 0;
	
#	include "plot_star.h"	
	
    }
}
