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
 * 6th order Adam-Bashford predictor with a 7th order Adam-Moulton corrector, 
 */


static point_2d	*tloc, *tv, *ta;

static int 	num_rk4;

static int	first_time = 1;

void restart_am7( int n, double *m, point_2d *x, point_2d *v )
{
    num_rk4 = 5;
}


void init_am7( int n, double *m, point_2d *x, point_2d *v )
{

    int		i;

    static int	last_n = -1;
    

    init_rk4( n, m, x, v );
    
    restart_am7( n, m, x, v );


    if( n != last_n )
    {
	last_n = n;
	
	/* allocate storage */
	if( !first_time )
	{
	    free( tloc );
	    free( tv );
	    free( ta );
	}
	
	tloc = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	tv = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );
	ta = (point_2d *)Safe_Malloc( sizeof( point_2d ) * n );

	for( i = 0; i < n; i++ )
	{
	    tloc[i].x = tloc[i].y = 0;
	    tv[i].x = tv[i].y = 0;
	    ta[i].x = ta[i].y = 0;
	}
    }

    first_time = 0;
	
}



static void collide_am7( int i, int j,
			point_2d *cur, point_2d *prev, point_2d *tloc,
			double *m,
			point_2d *vk[], point_2d *ak[],
			point_2d *tv, point_2d *ta,
			sys_param_type *sp, star_disp_type *sd
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
    tv[i].x = (m[i]*tv[i].x + m[j]*tv[j].x)*tmass;
    tv[i].y = (m[i]*tv[i].y + m[j]*tv[j].y)*tmass;

    
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
	tloc[i].x = (m[i]*tloc[i].x + m[j]*tloc[j].x)*tmass;
	tloc[i].y = (m[i]*tloc[i].y + m[j]*tloc[j].y)*tmass;
	
	
	for( k = 0; k < K; k++ )
	{
	    ak[k][i].x = (m[i]*ak[k][i].x + m[j]*ak[k][j].x)*tmass;
	    ak[k][i].y = (m[i]*ak[k][i].y + m[j]*ak[k][j].y)*tmass;
	}
	ta[i].x = (m[i]*ta[i].x + m[j]*ta[j].x)*tmass;
	ta[i].y = (m[i]*ta[i].y + m[j]*ta[j].y)*tmass;
	
	m[i] += m[j];
	if( m[i] >= collapsar )
	{
	    /* This _should_, by construction, never happen. */
	    if( verbose > 0 )
		fprintf( stderr, "%s:%d  combined mass too large: m[%d]=%g\n",
			__FILE__, __LINE__, i, m[i]/fm );
	    
	    m[i] = collapsar * .999;
	}
    }
    
    --sp->live_stars;
    sd->tstamp = time(NULL);
    set_buffer_factor( sd, sp );
    if( verbose > 1 )
	fprintf( stderr, "%s:%d  live_stars=%d  stars: %d,%d  m: %g+%g=%g\n", __FILE__, __LINE__, sp->live_stars, i, j, orig_m_i/fm, m[j]/fm, m[i]/fm );
    
    m[j] = 0.0;
    cur[j].x = cur[j].y = prev[j].x = prev[j].y = tloc[j].x = tloc[j].y = -1.0;
    for( k = 0; k < K; k++ )
    {
	vk[k][j].x = vk[k][j].y = 0;
	ak[k][j].x = ak[k][j].y = 0;
    }	    
    tv[j].x = tv[j].y = 0;
    ta[j].x = ta[j].y = 0;

    init_fna( max_stars, m, tloc, tv );
}


void move_stars_am7( point_2d *prev, point_2d *cur, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    double		m_i, m_j;
    double		previ_x, previ_y;
    double		tloci_x, tloci_y;
    double		ai_x, ai_y;

    register int        i, j, k;	/* star index */

    int			xpt, ypt;
    
    point_2d	*v_0 = vk[0];
    point_2d	*v_1 = vk[1];
    point_2d	*v_2 = vk[2];
    point_2d	*v_3 = vk[3];
    point_2d	*v_4 = vk[4];
    point_2d	*v_5 = vk[5];
    point_2d	*v_6 = vk[6];
    point_2d	*a_0 = ak[0];
    point_2d	*a_1 = ak[1];
    point_2d	*a_2 = ak[2];
    point_2d	*a_3 = ak[3];
    point_2d	*a_4 = ak[4];
    point_2d	*a_5 = ak[5];


    /*
     * we must use some other, self starting method for the first 7
     * calls to the Adam-Bashford predictor, Adam-Moulton corrector method
     */
    if( num_rk4 )
    {
	--num_rk4;
	move_stars_rk4( prev, cur, m, vk, ak, max_stars, sp, sd );
	return;
    }


    /*
     * predictor step
     */
    MARK( main_loop_p );
    for( i = 0; i < max_stars; i++ )
    {
	MARK( main_loop2_p );
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;
	
	previ_x = prev[i].x;
	previ_y = prev[i].y;
	ai_x = a_0[i].x;
	ai_y = a_0[i].y;
	
	MARK( force_calc_p );
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
		collide( i, j, cur, prev, m, vk, ak, sp, sd );
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
	
	/*
	 * Move the star
	 */
	
	MARK( move_star_p );
	if( m_i >= collapsar )
	{
	    /* I think it looks better if the collapsars don't move ;-) */
	    v_1[i].x = v_1[i].y = 0.0;
	    
	    a_0[i].x = a_0[i].y = 0.0;
	    ai_x = ai_y = 0.0;
	    a_1[i].x = a_1[i].y = 0.0;

	    tloc[i].x = prev[i].x;
	    tloc[i].y = prev[i].y;
	}
	else
	{
	    a_0[i].x = ai_x;
	    a_0[i].y = ai_y;

	    tv[i].x = v_1[i].x + (4227/1440.)*a_0[i].x - (7673/1440.)*a_1[i].x + (9482/1440.)*a_2[i].x - (6798/1440.)*a_3[i].x + (2627/1440.)*a_4[i].x - (425/1440.)*a_5[i].x;
	    tv[i].y = v_1[i].y + (4227/1440.)*a_0[i].y - (7673/1440.)*a_1[i].y + (9482/1440.)*a_2[i].y - (6798/1440.)*a_3[i].y + (2627/1440.)*a_4[i].y - (425/1440.)*a_5[i].y;
	    
	    tloc[i].x = prev[i].x + (4227/1440.)*v_1[i].x - (7673/1440.)*v_2[i].x + (9482/1440.)*v_3[i].x - (6798/1440.)*v_4[i].x + (2627/1440.)*v_5[i].x - (425/1440.)*v_6[i].x;
	    tloc[i].y = prev[i].y + (4227/1440.)*v_1[i].y - (7673/1440.)*v_2[i].y + (9482/1440.)*v_3[i].y - (6798/1440.)*v_4[i].y + (2627/1440.)*v_5[i].y - (425/1440.)*v_6[i].y;
	    
	}
    }


    /*
     * corrector step
     */
    MARK( main_loop_c );
    for( i = 0; i < max_stars; i++ )
    {
	MARK( main_loop2_c );
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;
	
	tloci_x = tloc[i].x;
	tloci_y = tloc[i].y;
	ai_x = ta[i].x;
	ai_y = ta[i].y;
	
	MARK( force_calc_c );
	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2, f;
	    register double f_d;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = tloc[j].x - tloci_x;
	    dy = tloc[j].y - tloci_y;
	    n2 = dx*dx + dy*dy;
	    
	    if( n2 < collide_dist )
	    {
		ta[i].x = ai_x;
		ta[i].y = ai_y;
		collide_am7( i, j, cur, prev, tloc, m, vk, ak, tv, ta, sp, sd );
		m_i = m[i];
		m_j = m[j];
		tloci_x = tloc[i].x;
		tloci_y = tloc[i].y;
		ai_x = ta[i].x;
		ai_y = ta[i].y;
	    }
	    else
	    {
		f = G / (sqrt(n2) * n2);
		
		f_d = dx * f;		
		ai_x += f_d * m_j;
		ta[j].x -= f_d * m_i;
		
		f_d = dy * f;
		ai_y += f_d * m_j;
		ta[j].y -= f_d * m_i;
	    }
	}
	
	/*
	 * Move the star
	 */
	
	MARK( move_star_c );
	if( m_i >= collapsar )
	{
	    /* I think it looks better if the collapsars don't move ;-) */
	    v_1[i].x = v_1[i].y = 0.0;
	    
	    ta[i].x = ta[i].y = 0.0;
	    ai_x = ai_y = 0.0;
	    a_1[i].x = a_1[i].y = 0.0;
	    
	    xpt = cur[i].x = tloc[i].x;
	    ypt = cur[i].y = tloc[i].y;
	}
	else
	{
	    ta[i].x = ai_x;
	    ta[i].y = ai_y;

	    v_0[i].x = v_1[i].x + (19087/60480.)*ta[i].x + (65112/60480.)*a_0[i].x - (46461/60480.)*a_1[i].x + (37504/60480.)*a_2[i].x - (20211/60480.)*a_3[i].x + (6312/60480.)*a_4[i].x - (863/60480.)*a_5[i].x;
	    v_0[i].y = v_1[i].y + (19087/60480.)*ta[i].y + (65112/60480.)*a_0[i].y - (46461/60480.)*a_1[i].y + (37504/60480.)*a_2[i].y - (20211/60480.)*a_3[i].y + (6312/60480.)*a_4[i].y - (863/60480.)*a_5[i].y;

	    xpt = cur[i].x = prev[i].x + (19087/60480.)*tv[i].x + (65112/60480.)*v_1[i].x - (46461/60480.)*v_2[i].x + (37504/60480.)*v_3[i].x - (20211/60480.)*v_4[i].x + (6312/60480.)*v_5[i].x - (863/60480.)*v_6[i].x;
	    ypt = cur[i].y = prev[i].y + (19087/60480.)*tv[i].y + (65112/60480.)*v_1[i].y - (46461/60480.)*v_2[i].y + (37504/60480.)*v_3[i].y - (20211/60480.)*v_4[i].y + (6312/60480.)*v_5[i].y - (863/60480.)*v_6[i].y;

	    vk[K-1][i].x = vk[K-1][i].y = 0;
	}
	ta[i].x = ta[i].y = 0;
	ak[K-1][i].x = ak[K-1][i].y = 0;
	
	MARK( plot_star );
#	include "plot_star.h"	
	
    }
}
