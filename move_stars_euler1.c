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
 * euler's method
 */


void restart_euler1( int n, double *m, point_2d *x, point_2d *v )
{
    return;
}



void init_euler1( int n, double *m, point_2d *x, point_2d *v )
{
    return;
}



void move_stars_euler1( point_2d *prev, point_2d *cur, double *m, 
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
    point_2d	*a_0 = ak[0];
    point_2d	*a_1 = ak[1];



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
	    v_0[i].x = v_1[i].x + ai_x;
	    xpt = cur[i].x = previ_x + v_1[i].x;
	    
	    
	    a_0[i].y = ai_y;
	    v_0[i].y = v_1[i].y + ai_y;
	    ypt = cur[i].y = previ_y + v_1[i].y;
	}
	ak[K-1][i].x = ak[K-1][i].y = 0;
	
	MARK( plot_star );
#	include "plot_star.h"	
	
    }
}
