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
 * Adam-Bashford method of 4th order
 */

static int 	num_rk4;

void restart_ab4( int n, double *m, point_2d *x, point_2d *v )
{
    num_rk4 = 4;
}



void init_ab4( int n, double *m, point_2d *x, point_2d *v )
{
    init_rk4( n, m, x, v );
    
    restart_ab4( n, m, x, v );
}



void move_stars_ab4( point_2d *prev, point_2d *cur, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    double		m_i, m_j;
    double		previ_x, previ_y;
    double		ai_x, ai_y;

    register int        i, j, k;	/* star index */

    int			xpt, ypt;
    
    point_2d	*v_0 = vk[0];
    point_2d	*v_1 = vk[1];
    point_2d	*v_2 = vk[2];
    point_2d	*v_3 = vk[3];
    point_2d	*v_4 = vk[4];
    point_2d	*a_0 = ak[0];
    point_2d	*a_1 = ak[1];
    point_2d	*a_2 = ak[2];
    point_2d	*a_3 = ak[3];


    /*
     * we must use some other, self starting method for the first 4
     * calls to the Adam-Bashford 4th order method
     */
    if( num_rk4 )
    {
	--num_rk4;
	move_stars_rk4( prev, cur, m, vk, ak, max_stars, sp, sd );
	return;
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
	
	MARK( move_star );
	if( m_i >= collapsar )
	{
	    /* I think it looks better if the collapsars don't move ;-) */
	    v_1[i].x = v_1[i].y = 0.0;
	    
	    a_0[i].x = a_0[i].y = 0.0;
	    ai_x = ai_y = 0.0;
	    a_1[i].x = a_1[i].y = 0.0;
	    
	    xpt = cur[i].x = prev[i].x;
	    ypt = cur[i].y = prev[i].y;
	}
	else
	{
	    a_0[i].x = ai_x;
	    a_0[i].y = ai_y;

	    v_0[i].x = v_1[i].x + (55/24.)*a_0[i].x - (59/24.)*a_1[i].x + (37/24.)*a_2[i].x - (9/24.)*a_3[i].x;
	    v_0[i].y = v_1[i].y + (55/24.)*a_0[i].y - (59/24.)*a_1[i].y + (37/24.)*a_2[i].y - (9/24.)*a_3[i].y;
	    
            xpt = cur[i].x = previ_x + (55/24.)*v_1[i].x - (59/24.)*v_2[i].x + (37/24.)*v_3[i].x - (9/24.)*v_4[i].x;
            ypt = cur[i].y = previ_y + (55/24.)*v_1[i].y - (59/24.)*v_2[i].y + (37/24.)*v_3[i].y - (9/24.)*v_4[i].y;


	    vk[K-1][i].x = vk[K-1][i].y = 0;
	}
	ak[K-1][i].x = ak[K-1][i].y = 0;
	
	MARK( plot_star );
#	include "plot_star.h"	
	
    }
}
