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
 * this file contains a bunch of routines that deal with the constants
 * of motion that all star systems have.  (total energy, angular momentum,
 * linear momentum, etc)
 */


double kinetic_energy( int i, int n, point_2d *cur, double *m, point_2d *v )
{
    return .5 * m[i] * (v[i].x*v[i].x + v[i].y*v[i].y);
}


double binding_energy( int i, int n, point_2d *cur, double *m, point_2d *v )
{
    double	dx, dy;
    double	sum;
    int		j;

    sum = 0;
    for( j = 0; j < n; j++ )
    {
	if( m[j] <= 0. )
	    continue;

	if( i == j )
	    continue;

	dx = cur[i].x - cur[j].x;
	dy = cur[i].y - cur[j].y;
	sum += m[j]/sqrt( dx*dx + dy*dy );
    }

    return -.5 * G * m[i] * sum;
}

    


void calc_sys_const( sys_const *sc,
		    int n, point_2d *cur, double *m,  point_2d *v )
{
    int		i;
    double	m_i;
    double	k, u;
    
    sc->k = 0;
    sc->u = 0;
    sc->e = 0;
    sc->cm.x = 0;
    sc->cm.y = 0;
    sc->p.x = 0;
    sc->p.y = 0;
    sc->l = 0;
    sc->m = 0;

    for( i = 0; i < n; i++ )
    {
	m_i = m[i];
	if( m_i <= 0. )
	    continue;

	k = kinetic_energy( i, n, cur, m, v );
	u = binding_energy( i, n, cur, m, v );

	sc->k += k;

	sc->u += u;
	sc->e += k + u;

	sc->cm.x += m_i * cur[i].x;
	sc->cm.y += m_i * cur[i].y;

	sc->p.x += m_i * v[i].x;
	sc->p.y += m_i * v[i].y;

	sc->l += m_i * (cur[i].x * v[i].y - cur[i].y * v[i].x);

	sc->m += m_i;
    }

    sc->cm.x /= sc->m;
    sc->cm.y /= sc->m;

}


void do_prt_sys_const( int n, point_2d *cur, double *m,  point_2d *v, char *file, int line )
{
    sys_const	sc;

    calc_sys_const( &sc, n, cur, m, v );

    fprintf( stderr, "%s:%d  k=%9.3e + u=%9.3e = %9.3e  l=%8.5g  m=%7.4g\n", file, line, sc.k*fv*fv/fm, sc.u*fv*fv/fm, (sc.k + sc.u)*fv*fv/fm, sc.l*fv/fm, sc.m/fm );
    fprintf( stderr, "%s:%d  cm=(%6g,%6g)  p=(%6g,%6g)\n", file, line, sc.cm.x, sc.cm.y, sc.p.x*fv/fm, sc.p.y*fv/fm );
}



#define per_diff(a,b) ( (a) ? 100.*((a)-(b))/(a) : (a)-(b) )

void do_prt_sys_const_delta( sys_const sc0, int n, point_2d *cur, double *m,  point_2d *v, char *file, int line )
{
    sys_const	sc1;

    calc_sys_const( &sc1, n, cur, m, v );

    fprintf( stderr, "%s:%d  k + u = %9.3e%%  l=%8.5g%%  m=%8.5g%%\n", file, line, per_diff(sc0.k + sc0.u,sc1.k + sc1.u), per_diff(sc0.l,sc1.l), per_diff(sc0.m,sc1.m) );
    fprintf( stderr, "%s:%d  cm=(%6g,%6g)  p=(%6g,%6g)\n", file, line, sc0.cm.x-sc1.cm.x, sc0.cm.y-sc1.cm.y, sc0.p.x-sc1.p.x, sc0.p.y-sc1.p.y );
}

