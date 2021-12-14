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



void init_sys_4( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp )
{
    int		i, k;

    sp->star_distrib     = 0;
    sp->size             = 2.1400000000000000e+02;
    sp->num_collapsar    = 0;
    sp->mono_stars       = 4;

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

    sp->min_stars        = 1;
    sp->num_add          = 0;
    sp->live_stars       = 4;

/* star 0 */

    m[0]        = 2.7865168539325840e+00;
    cur[0].x    = 1.5835647927008775e+01;
    cur[0].y    = 1.8923118463283608e+01;
    prev[0].x   = 1.5835647927008775e+01;
    prev[0].y   = 1.8923118463283608e+01;
    vk[0][0].x  = 9.0941882306722933e-03;
    vk[0][0].y  = 1.6257009059518187e-03;
    ak[0][0].x  = 6.9931581889377990e-06;
    ak[0][0].y  = -8.0395052303351502e-06;

/* star 1 */

    m[1]        = 5.4831460674157304e+00;
    cur[1].x    = 6.0521453125534116e+01;
    cur[1].y    = 9.4644832984345086e+00;
    prev[1].x   = 6.0521453125534116e+01;
    prev[1].y   = 9.4644832984345086e+00;
    vk[0][1].x  = 2.5826514208378845e-03;
    vk[0][1].y  = -1.4348317268012971e-02;
    ak[0][1].x  = -8.3617935918065052e-06;
    ak[0][1].y  = -4.3824681501481231e-07;

/* star 2 */

    m[2]        = 2.5168539325842696e+00;
    cur[2].x    = 8.3340088264855723e+00;
    cur[2].y    = -2.3248024535545774e+01;
    prev[2].x   = 8.3340088264855723e+00;
    prev[2].y   = -2.3248024535545774e+01;
    vk[0][2].x  = -1.4385862343579154e-02;
    vk[0][2].y  = -1.3381252605483343e-03;
    ak[0][2].x  = 3.0889420909050185e-06;
    ak[0][2].y  = 9.3058808330760419e-06;

/* star 3 */

    m[3]        = 5.2134831460674160e+00;
    cur[3].x    = -7.6139068509249114e+01;
    cur[3].y    = -8.8449218374312917e+00;
    prev[3].x   = -7.6139068509249114e+01;
    prev[3].y   = -8.8449218374312917e+00;
    vk[0][3].x  = -6.3202458968509645e-04;
    vk[0][3].y  = 1.4867554009648933e-02;
    ak[0][3].x  = 3.5653642534100792e-06;
    ak[0][3].y  = 2.6539749189938063e-07;


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
