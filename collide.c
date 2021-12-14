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
 * This routine is called (via the macro collide()) when two stars get
 * close enough that they should merge.  The minimal separation
 * distance is important.  If you don't try to collapse two nearby
 * stars together soon enough, they reach very high spin velocities.
 * Then, the slightest disturbance can cause them to step over and
 * past each other in a single movement.  In reality they would have
 * had to have passed through each other to do this last step so they
 * should have collided.  After this last step, they are now far
 * enough apart and have such high velocities, that their attraction
 * can no longer keep them together and they fly off at very high
 * speeds.
 *
 * The smaller the FV value, the larger this value needs to be,
 * because movement is coarser and each step can be larger.  The collide
 * distance can be changed with the -l command line option.
 *
 * Be aware that besides this overstep artifact caused by discrete
 * movement, there is also the sling shot effect which is quite real.
 * The sling shot effect is actually used for interplanetary probes
 * and such.  It works by transferring some of the energy from a
 * larger body (planet, star) to the smaller body (probe), allowing it
 * to gain a fair amount of speed.  You can usually tell the
 * difference between the results of the two phenomenon, because the
 * overstep artifact shoots stars off at a very high speed, and the
 * sling shot effect only shoots stars off at a high speed.
 *
 * Now that I can calculate the star systems constants of motion, you
 * can always distinguish these two phenomenons by checking the total
 * energy of the system.  If the energy has increased, then it is the
 * overstep artifact.  If the energy has remained constant, then it is
 * the sling shot effect.
 *
 *
 * Since collapsars don't move, you will never see the sling shot
 * effect from a collapsar.  Because collapsars are heavy, you will
 * often see the overstep artifact from them.
 */

void do_collide( int i, int j, point_2d *cur, point_2d *prev, double *m,
		point_2d *vk[], point_2d *ak[],
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
	fprintf( stderr, "%s:%d  live_stars=%d  stars: %d,%d  m: %g+%g=%g\n", file, line, sp->live_stars, i, j, orig_m_i/fm, m[j]/fm, m[i]/fm );
    
    m[j] = 0.0;
    cur[j].x = cur[j].y = prev[j].x = prev[j].y = -1.0;
    for( k = 0; k < K; k++ )
    {
	vk[k][j].x = vk[k][j].y = 0;
	ak[k][j].x = ak[k][j].y = 0;
    }	    

    init_fna( max_stars, m, prev, vk[1] );
}


