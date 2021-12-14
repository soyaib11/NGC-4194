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



void init_sys_1c( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp )
{
    int		i, k;

    sp->star_distrib     = 0;
    sp->size             = 2.1400000000000000e+02;
    sp->num_collapsar    = 2;
    sp->mono_stars       = 1;

    sp->star_circle      = 0;
    sp->cir_dist         = 0;

    sp->star_line        = 0;
    sp->rnd_spacing      = 0;

    sp->do_bounce        = 0;
    sp->no_speed         = 0;
    sp->few_stars        = 0;
    sp->drift            = 0;
    sp->min_angular_mom  = 0;

    sp->energy_fact      = 0;
    sp->calc_energy_fact = 0;

    sp->min_stars        = 0;
    sp->num_add          = 0;
    sp->live_stars       = 1;

/* star 0 */

    m[0]        = 4.0000000000000000e+00;
    cur[0].x    = 7.0000000000000000e+01;
    cur[0].y    = 0.0000000000000000e+00;
    prev[0].x   = 7.0000000000000000e+01;
    prev[0].y   = 0.0000000000000000e+00;
    vk[0][0].x  = 0.0000000000000000e+00;
    vk[0][0].y  = 1.5000000000000000e-02;
    ak[0][0].x  = 6.9931581889377999e-06;
    ak[0][0].y  = -8.0395052303351518e-06;

/* star 1 */

    m[1]        = 1.6001600000000000e+01;
    cur[1].x    = 0.0000000000000000e+00;
    cur[1].y    = 0.0000000000000000e+00;
    prev[1].x   = 0.0000000000000000e+00;
    prev[1].y   = 0.0000000000000000e+00;


/* star 2 */

    m[2]        = 1.6001600000000000e+01;
    cur[2].x    = 3.0000000000000000e+00;
    cur[2].y    = 0.0000000000000000e+00;
    prev[2].x   = 0.0000000000000000e+00;
    prev[2].y   = 0.0000000000000000e+00;


    /* scale things to the current accuracy */

    for( i = 0; i < max_stars; i++ )
    {
	if( m[i] <= 0 )
	    continue;
	
	m[i] *= fm;

	for( k = 0; k < K; k++ )
	{
	    vk[k][i].x *= fv_inv;
	    vk[k][i].y *= fv_inv;

	    ak[k][i].x *= fv_inv*fv_inv;
	    ak[k][i].y *= fv_inv*fv_inv;
	}
    }
}
