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



#include <time.h>

/*
 * Main driver routine for initializing the star system.
 */

void init_stars( point_2d *cur, point_2d *prev, double *m,
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    int		i, j, k;

    point_2d	*v_0 = vk[0];

    int		tmp_color;


    /* clean up from last time */
    sd->tstamp = time(NULL);
    if( verbose > 1 || (verbose == 1 && last_init == 0) )
    {
	fprintf( stderr, "                 started with                |             ended with\n" );
	fprintf( stderr, "            gravity                     no   | life\n" );
	fprintf( stderr, "stars/bin/++  well  size  cir lin bnce speed | span  collide  edge  ++    steps\n" );
    }

    if( verbose > 0 && last_init != 0 )
    {
	fprintf( stderr,
		"%5d/%2d/%2d %4d %8.2f %2d %3d %3d %5d   | %4ld %5d %7d %4d%9d\n",
		orig_sp.live_stars, orig_sp.live_stars-sp->mono_stars,
		orig_sp.num_add, sp->num_collapsar,
		sp->size, sp->star_circle, sp->star_line, sp->do_bounce,
		sp->no_speed,
		sd->tstamp - last_init, num_collide, num_edge,
		sp->num_add, sd->num_steps );

	if( verbose > 3 )
	{
	    prt_sys_const( max_stars, cur, m, v_0 );
	    prt_sys_const_delta( orig_sc, max_stars, cur, m, v_0 );
	}
    }

    if( last_init != 0 )
    {
	update_screen( sp, sd, &display );
	XFlush( display.dpy );
    }
    clear_disp_pt( sd );

    last_init = sd->tstamp;
    num_collide = 0;
    num_edge = 0;


    /* make sure the array is clean at the end... */
    for( i = 0; i < max_stars; i++ )
    {
	prev[i].x = cur[i].x = 0;
	prev[i].y = cur[i].y = 0;
	m[i] = 0;
	sd->last_disp[i].x = far_dist;
	sd->last_disp[i].y = far_dist;

	for( k = 0; k < K; k++ )
	{
	    vk[k][i].x = vk[k][i].y = 0;
	    ak[k][i].x = ak[k][i].y = 0;
	}	    
    }


#if 1
    /* setup up the system parameters */
    set_sys_param( sp, max_stars );

    /* set the mass, locations, velocities and accelerations */
    set_xmva( cur, prev, m, vk, ak, max_stars, sp );
#else
    /* use a predefined star system */
/*    init_sys_1a( cur, prev, m, vk, ak, max_stars, sp ); */
/*    init_sys_4( cur, prev, m, vk, ak, max_stars, sp ); */
/*    sd->color_number = 0; */
    init_sys_8( cur, prev, m, vk, ak, max_stars, sp );
/*    init_sys_8b( cur, prev, m, vk, ak, max_stars, sp ); */
#endif
    
    /* remember the original system parameters */
    orig_sp = *sp;


    /* set up the star_disp variables */
    set_star_disp( sd, sp );
    set_buffer_factor( sd, sp );


    /* do any movement initialization that is required */
    calc_sys_const( &orig_sc, max_stars, cur, m, v_0 );
    if( (orig_sc.k + orig_sc.u)*fv*fv/fm > -1E-5 )
    {
	fprintf( stderr, "%s:%d  Error:  too much energy, system is not bound.\n", __FILE__, __LINE__ );
	prt_sys_const( max_stars, cur, m, vk[0] );
	prt_sys_const_delta( orig_sc, max_stars, cur, m, vk[0] );
    }
    init_fna( max_stars, m, cur, v_0 );


    /*
     * set up a random star color system for multi-color mode
     */
    for( i = 0; i < NUM_COLORS; i++ )
	sd->rnd_colors[ i ] = i;
    
    for( i = 0; i < NUM_COLORS-1; i++ )
    {
	j = i + lrand48() % (NUM_COLORS - i);
	tmp_color = sd->rnd_colors[ i ];
	sd->rnd_colors[ i ] = sd->rnd_colors[ j ];
	sd->rnd_colors[ j ] = tmp_color;
    }

    for( i = 0; i < max_stars; i++ )
	sd->star_color[ i ] = sd->rnd_colors[ i % NUM_COLORS ];
    


    /*
     * print out the star configuration
     */

    if( verbose > 2 )
    {
	double	u, k;
	double	dx, dy;
	double	tdist;
	
	fprintf( stderr, "star   location         velocity     mass   tdist  kinetic   binding   -k/u\n" );
    
	for( i = 0; i < sp->live_stars; i++ )
	{
	    dx = cur[i].x;
	    dy = cur[i].y;
	    tdist = sqrt( dx*dx + dy*dy );
	    
	    k = kinetic_energy( i, max_stars, cur, m, v_0 );
	    u = binding_energy( i, max_stars, cur, m, v_0 );
	    
	    fprintf( stderr,
		    "%2d  %6.1f,%6.1f  %7.4f,%7.4f %6.3f  %6.2f %9.6f %9.6f %9.6f\n",
		    i, cur[i].x, cur[i].y, fv*v_0[i].x, fv*v_0[i].y, m[i]/fm,
		    tdist, k*fv*fv/fm, u*fv*fv/fm, -k/u );
	}
    }
#if 0
    /* dump the star system configuration for later reuse */
    dump_sys( cur, prev, m, vk, ak, max_stars, sp );
#endif
    

    /* clean up the display */
    sleep(1);
    XFillRectangle(display.dpy, display.win, display.erase_gc,
		   0,0, winW, winH);
    plot_collapsars( cur, m, sp, sd, &display );
    XFlush( display.dpy );
}



/*
 * add a new star to the system or reinitialize it
 */

/* note:  set_xmva() has similar functionality... change both places */

void new_star( point_2d *cur, point_2d *prev, double *m, 
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp, star_disp_type *sd,
	      int add_collapsar )
{
    int		i, j;
    int		found;
    double	dx, dy;
    double	tdist;

    double	dist, cur_dist, max_dist;
    double	theta;
    


    static int	check = 0;


    /*
     * find a good slot to use.  (rotate through so we get different colors)
     */
    for( i = 0, found = -1; i < max_stars; i++ )
    {
	check++;
	if( check >= max_stars )
	    check = check % max_stars;
	
	if( m[check] <= 0 )
	{
	    found = check;
	    break;
	}
    }
    
    if( (sp->num_add)-- <= 0 || found == -1 )
    {
	init_stars( prev, cur, m, vk, ak, max_stars, sp, sd );
	return;
    }
    

    /*
     * create the star
     */

    max_dist = far_dist * far_dist;
    do
    {
	/*
	 * use appropriate distribution of stars
	 */
	theta = 2*M_PI * drand48();
	switch( sp->star_distrib )
	{
	case DIST_CENTER:
	    dist = drand48() * sp->size;
	    break;
	    
	case DIST_UNIFORM:
	    dist = sqrt( drand48() ) * sp->size;
	    break;
	    
	case DIST_RING:
	    dist = sqrt(sqrt( drand48() )) * sp->size;
	    break;
	    
	default:
	    /* this shouldn't happen... */
	    dist = 100;
	}

	if( add_collapsar && !sp->do_bounce )
	{
	    prev[found].x = cur[found].x = dist * .5 * cos( theta );
	    prev[found].y = cur[found].y = dist * .5 * sin( theta );
	}
	else
	{
	    prev[found].x = cur[found].x = dist * cos( theta );
	    prev[found].y = cur[found].y = dist * sin( theta );
	}
	
	
	/*
	 * don't put stars too close together
	 */
	cur_dist = far_dist * far_dist;
	for( j = 0; j < max_stars; j++ )
	{
	    if( found == j )
		continue;

	    if( m[j] <= 0 )
		continue;
	    
	    dx = prev[found].x - prev[j].x;
	    dy = prev[found].y - prev[j].y;
	    dist = dx*dx + dy*dy;
	    
	    if( dist < cur_dist )
		cur_dist = dist;
	}
    }
    while( cur_dist < 30*30
	  && ( cur_dist < 20*20 || drand48() < .5 )
	  && ( cur_dist < 10*10 || drand48() < .75 )
	  );
    
    if( add_collapsar )
    {
	m[found] = collapsar * 1.0001;
	vk[0][found].x = vk[0][found].y = ak[0][found].x = ak[0][found].y = 0;
    }
    else
    {
	if( sp->num_collapsar )
	    m[found] = fm*.1;
	else
	    m[found] = fm*.5;

	set_va( cur, prev, m, vk, ak, sp, found );
    }

    sd->last_disp[found].x = far_dist;
    sd->last_disp[found].y = far_dist;


    (sp->live_stars)++;
    (sd->live_erase)++;
    sd->tstamp = time(NULL);
    set_buffer_factor( sd, sp );

    if( verbose > 2 )
    {
	dx = cur[found].x;
	dy = cur[found].y;
	tdist = sqrt( dx*dx + dy*dy );
	
	fprintf( stderr,
		"%s:%d  star %2d @ (%6.1f,%6.1f)=|%7.4f,%7.4f|*%6.3f  tdist=%6.2f\n",
		__FILE__, __LINE__, found, cur[found].x, cur[found].y, fv*vk[0][found].x,
		fv*vk[0][found].y, m[found]/fm, tdist );
    }


    /* do any movement initialization that is required */
    init_fna( max_stars, m, cur, vk[0] );

    /* make sure collapsars get plotted... */
    if( add_collapsar )
	plot_collapsars( cur, m, sp, sd, &display );
}

