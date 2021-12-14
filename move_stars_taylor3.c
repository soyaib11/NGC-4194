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
 * taylor series expansion to the 2nd order with the 3rd order approximated
 * by using the derivative the the 2nd degree Lagrange polynomial.
 */


static int 	num_rk4;

void restart_taylor3( int n, double *m, point_2d *x, point_2d *v )
{
    num_rk4 = 3;
}



void init_taylor3( int n, double *m, point_2d *x, point_2d *v )
{
    init_rk4( n, m, x, v );
    
    restart_taylor3( n, m, x, v );
}



void move_stars_taylor3( point_2d *prev, point_2d *cur, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    double		m_i, m_j;
    double		previ_x, previ_y;
    double		ai_x, ai_y;
    double		Dt_a;		/* derivative of the acceleration */

    register int        i, j, k;	/* star index */

    int			xpt, ypt;
    
    point_2d	*v_0 = vk[0];
    point_2d	*v_1 = vk[1];
    point_2d	*a_0 = ak[0];
    point_2d	*a_1 = ak[1];
    point_2d	*a_2 = ak[2];



    /*
     * we must use some other, self starting method for the first 3
     * calls to the taylor method
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
	 *
	 * You want each step to be as small as possible and this is done
	 * by making FV as large as possible.  If the steps are small
	 * enough, then we are approximating infinitesimal movements.
	 *
	 * The movement is based off of the formula:  (with t=1)
	 * x(t) = x0 + v t + 1/2 a t^2 + 1/(2*3) (Dt a) t^3
	 * v(t) = v0 + a t + 1/2 (Dt a) t^2
	 *
	 * In theory, these formulas should include an infinite number
	 * of derivatives of x(t), but the additional terms are not
	 * easy to calculate.  Even Dt a is just an estimate.
	 *
	 *
	 * note:
	 * I am currently using the derivative of the 2nd degree
	 * Lagrange polynomial to estimate Dt_a.  I tried using the
	 * 3rd degree Lagrange Polynomial, but it didn't work as well.
	 * Polynomials of high degree tend to not be very smooth so
	 * their derivatives can make things worse.
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
/*	    Dt_a = (11./12)*ai_x - 1.5*a_1[i].x + .75*a_2[i].x
		- (1./6)*a_3[i].x; */
	    Dt_a = (3./4)*ai_x - a_1[i].x + (1./4)*a_2[i].x;
	    v_0[i].x = v_1[i].x + ai_x + Dt_a;
	    xpt = cur[i].x = previ_x + (v_1[i].x + ((1./2)*ai_x
						    + (1./3)*Dt_a));
	    
	    
	    a_0[i].y = ai_y;
/*	    Dt_a = (11./12)*ai_y - 1.5*a_1[i].y + .75*a_2[i].y
		- (1./6)*a_3[i].y; */
	    Dt_a = (3./4)*ai_y - a_1[i].y + (1./4)*a_2[i].y;
	    v_0[i].y = v_1[i].y + ai_y + Dt_a;
	    ypt = cur[i].y = previ_y + (v_1[i].y + ((1./2)*ai_y
						    + (1./3)*Dt_a));
	}
	ak[K-1][i].x = ak[K-1][i].y = 0;
	
	MARK( plot_star );
#	include "plot_star.h"	
	
    }
}
