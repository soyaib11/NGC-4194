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
 * This routine dumps a star system into a format suitable for including
 * in another C program.  This routine was used to create the init_sys_*()
 * routines.
 */

void dump_sys( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp )
{
    int		i, k;
    
    /*
     * dump the variables so that we can recreate this system later.
     */

    fprintf( stderr, "    sp->star_distrib     = %d;\n", sp->star_distrib );
    fprintf( stderr, "    sp->size             = %.16e;\n", sp->size );
    fprintf( stderr, "    sp->num_collapsar    = %d;\n", sp->num_collapsar );
    fprintf( stderr, "    sp->mono_stars       = %d;\n", sp->mono_stars );
    fprintf( stderr, "\n" );
    fprintf( stderr, "    sp->star_circle      = %d;\n", sp->star_circle );
    fprintf( stderr, "    sp->cir_dist         = %.16e;\n", sp->cir_dist );
    fprintf( stderr, "\n" );
    fprintf( stderr, "    sp->star_line        = %d;\n", sp->star_line );
    fprintf( stderr, "    sp->rnd_spacing      = %d;\n", sp->rnd_spacing );
    fprintf( stderr, "\n" );
    fprintf( stderr, "    sp->do_bounce        = %d;\n", sp->do_bounce );
    fprintf( stderr, "    sp->no_speed         = %d;\n", sp->no_speed );
    fprintf( stderr, "    sp->few_stars        = %d;\n", sp->few_stars );
    fprintf( stderr, "    sp->drift            = %d;\n", sp->drift );
    fprintf( stderr, "    sp->min_angular_mom  = %d;\n", sp->min_angular_mom );
    fprintf( stderr, "\n" );
    fprintf( stderr, "    sp->energy_fact      = %.16e;\n", sp->energy_fact );
    fprintf( stderr, "    sp->calc_energy_fact = %.16e;\n", sp->calc_energy_fact );
    fprintf( stderr, "\n" );
    fprintf( stderr, "    sp->min_stars        = %d;\n", sp->min_stars );
    fprintf( stderr, "    sp->num_add          = %d;\n", sp->num_add );
    fprintf( stderr, "    sp->live_stars       = %d;\n", sp->live_stars );
    
    for( i = 0; i < max_stars; i++ )
    {
	if( m[i] <= 0 )
	    continue;
	
	fprintf( stderr, "\n/* star %d */\n", i );
	fprintf( stderr, "\n    m[%d]        = %.16e;\n", i, m[i]/fm );
	fprintf( stderr, "    cur[%d].x    = %.16e;\n", i, cur[i].x );
	fprintf( stderr, "    cur[%d].y    = %.16e;\n", i, cur[i].y );
	fprintf( stderr, "    prev[%d].x   = %.16e;\n", i, prev[i].x );
	fprintf( stderr, "    prev[%d].y   = %.16e;\n", i, prev[i].y );

	for( k = 0; k < K; k++ )
	{
	    if( vk[k][i].x != 0 || vk[k][i].y != 0 )
	    {
		fprintf( stderr, "    vk[%d][%d].x  = %.16e;\n",
			k, i, vk[k][i].x*fv );
		fprintf( stderr, "    vk[%d][%d].y  = %.16e;\n",
			k, i, vk[k][i].y*fv );
	    }
	    
	    if( ak[k][i].x != 0 || ak[k][i].y != 0 )
	    {
		fprintf( stderr, "    ak[%d][%d].x  = %.16e;\n",
			k, i, ak[k][i].x*fv*fv );
		fprintf( stderr, "    ak[%d][%d].y  = %.16e;\n",
			k, i, ak[k][i].y*fv*fv );
	    }
	    
	}
    }
}

