#include "xstar.h"
#include "xstar_ext.h"
#include <sys/time.h>
#include <memory.h>


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com) and others
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/



/*
 * this is the main loop that animates the stars.  It controls the creation
 * of the star systems, the movement, and all user commands.
 */

void
Animate()
{
    int			num_erase;
    int			bind_check;
    point_2d		*cur, *prev;	/* current and previous positions */
    point_2d		*temp_points;
    point_2d		*vk[K];		/* velocities		*/
    double              *m;             /* masses */
    point_2d		*ak[K];		/* accelerations	*/
    XEvent              xev;
    int			found;
    int			i, k;
    double		found_mass, tmass, dist;

    sys_param_type	sp;
    star_disp_type	sd;

    int			color_skip = 0;

    int			paused = FALSE;

    memset( &sp, 0, sizeof( sp ) );
    memset( &sd, 0, sizeof( sd ) );
    sd.next_free = 1;


    /* Allocate memory. */
    if( num_disp_pt )
    {
	sd.max_points = 512;		/* replotting screen needs large buf */

	sd.disp_pts = (disp_point_type *)Safe_Malloc( sizeof( disp_point_type )
						  * num_disp_pt );

	sd.hash_index = (ushort_t *)Safe_Malloc( sizeof( ushort_t )
					     * HASH_TABLE_SIZE );

    }
    else
    {
	sd.max_points = max_stars*10;	/* will rarely be more than 1/4 full */
	sd.tmp_pts = (XPoint *) Safe_Malloc(sizeof(XPoint) * sd.max_points);
	sd.pixels = (int *) Safe_Malloc(sizeof(int) * sd.max_points);
    }
    sd.points = (XPoint *) Safe_Malloc(sizeof(XPoint) * sd.max_points);
    cur = (point_2d *) Safe_Malloc( sizeof( point_2d ) * max_stars );
    prev = (point_2d *) Safe_Malloc( sizeof( point_2d ) * max_stars );
    sd.last_disp = (XPoint *) Safe_Malloc( sizeof( XPoint ) * max_stars );
    m = (double *) Safe_Malloc(sizeof(double) * max_stars);
    sd.star_color = (int *)Safe_Malloc( sizeof( int ) * max_stars );

    for( k = 0; k < K; k++ )
    {
	vk[k] = (point_2d *)Safe_Malloc( sizeof( point_2d ) * max_stars );
	ak[k] = (point_2d *)Safe_Malloc( sizeof( point_2d ) * max_stars );
    }


    /* set up colors */
    sd.color_number = lrand48() % NUM_COLORS;
    if( rotate_colors || multi_colors)
	init_colors( colors, color_gcs, &sd, NUM_COLORS );


    /* initialize the system */
    set_redraw_full( &sd );
    clear_disp_pt( &sd );
    bind_check = 0;
    num_erase = 0;
    
    init_stars( prev, cur, m, vk, ak, max_stars, &sp, &sd );


    /*
     * Seemingly endless loop.
     */
    while( TRUE )
    {
	sd.raw_num_steps++;
	sd.raw_hsteps++;
	sd.num_steps = sd.raw_num_steps * sd.hstep_scale;
	sd.hsteps = sd.raw_hsteps * sd.hstep_scale;
	if( sd.hsteps >= MAX_HSTEP - 128 )
	{
	    reset_hstep( &sd );
	    update_screen( &sp, &sd, &display );
	}

	sd.num_visible = 0;
	
	/* Age the the points. */
	temp_points = prev;
	prev = cur;
	cur = temp_points;
	
	temp_points = vk[K-1];
	for( k = K-1; k > 0; k-- )
	    vk[k] = vk[k-1];
	vk[0] = temp_points;
	
	temp_points = ak[K-1];
	for( k = K-1; k > 0; k-- )
	    ak[k] = ak[k-1];
	ak[0] = temp_points;
	

	/* calculate and display the new locations of each star */
	if( sp.do_bounce )
	    check_bounce( prev, cur, m, vk, ak, max_stars, &sp, &sd );

	move_fna( prev, cur, m, vk, ak, max_stars, &sp, &sd );

	
	/* 
	 * if there are not enough stars to be interesting, then add
	 * a new star or start over with a new system.
	 */
	if( sp.live_stars <= sp.min_stars
	   || ( sp.num_add > 1 && sp.live_stars <= sp.min_stars + 1 )
	   )
	{
	    new_star( prev, cur, m, vk, ak, max_stars, &sp, &sd, 0 );
	    continue;
	}
	

        /*
	 * Display stars
	 *
	 * some of the terms in these formulas should really be adjustable
	 * depending on the speed of the computer and the quality of the
	 * X server and compiler.  max_disp_skip in particular.
	 */
	sd.num_disp_skipped++;
	if( num_disp_pt )
	{
	    int		pts_to_plot = DIFF( sd.next_free, sd.next_disp );
	    
	    if( pts_to_plot > sd.num_visible
	       && ( sd.num_disp_skipped > sd.buffer_factor /* long time between disps */
		   || pts_to_plot >= 4*sd.num_visible		/* jerky lines */
		   || sp.live_stars >= DIFF( sd.next_erase, sd.next_free )	/* buf full */
		   )
	       )
		update_screen( &sp, &sd, &display );
	}
	else
	{
	    if( sd.points_used > sd.num_visible
	       && ( sd.num_disp_skipped > sd.buffer_factor /* long time between disps */
		   || sd.points_used >= 5*sd.num_visible	/* jerky lines */
		   || sd.points_used + stars + MAX_COLLAPSAR >= sd.max_points	/* buf full */
		   )
	       )
		update_screen( &sp, &sd, &display );
	}	
	
        /*
	 * Check for events.
	 */
	sd.num_poll_skipped++;
        if( sd.num_poll_skipped >= POLL_SKIP )
        {
	    int		check_event;
	    
	    sd.num_poll_skipped = 0;
	
	    check_event = XEventsQueued( display.dpy, QueuedAfterFlush );
	    if( !check_event )
	    {
		/* Delay so we don't use all of the cpu time. */
#ifdef USE_USLEEP
		if( delay != 0 ) usleep(delay);
		if( paused ) usleep(3000*1000);
#else
		if( delay != 0 ) nap(0,delay);
		if( paused ) nap(3,0);
#endif
	    }
	
	
	    /* check to see if the star system is still bound */
	    if( verbose > 1 && (bind_check++)*fv_inv > 128 )
	    {
		sys_const	sc;
		
		bind_check = 0;
		
		calc_sys_const( &sc, max_stars, cur, m, vk[0] );
		if( (sc.k + sc.u)*fv*fv/fm > -1E-7 )
		{
		    fprintf( stderr, "%s:%d  Error:  too much energy, system is not bound.\n", __FILE__, __LINE__ );
		    prt_sys_const( max_stars, cur, m, vk[0] );
		    prt_sys_const_delta( orig_sc, max_stars, cur, m, vk[0] );
		}
		else if( verbose > 3 )
		{
/*		    prt_sys_const( max_stars, cur, m, vk[0] ); */
		    prt_sys_const_delta( orig_sc, max_stars, cur, m, vk[0] );
		}
	    }


	    /* change the color of the stars */
 	    if( (color_skip++)*fv_inv > 80000./(POLL_SKIP*NUM_COLORS) )
	    {
		color_skip = 0;
		sd.color_number = (sd.color_number + 1) % NUM_COLORS;

		plot_collapsars( cur, m, &sp, &sd, &display );
	    }
	
	
	    /* Screen saver/background:  restart after 7 minutes
	     * with out a star dieing */
	    if( timeout || root )
	    {
		double	lost_stars;
		
		time_t	cur_time = time(NULL);
		
		if( cur_time - sd.tstamp > 420*( 1 + sp.do_bounce ) )
		    init_stars( prev, cur, m, vk, ak, max_stars, &sp, &sd );
		else if( cur_time - sd.tstamp > 210*( 1 + sp.do_bounce )
			&& sp.num_add > 1 )
		    new_star( prev, cur, m, vk, ak, max_stars, &sp, &sd, 0 );

		
		if( sd.live_erase < sp.live_stars )
		    sd.live_erase = sp.live_stars;
		
		/*
		 * should we clear the screen?
		 */
		if( !num_disp_pt && !sp.star_line )
		    if( sp.do_bounce )
		    {
			if( sd.erase_disps*fv_inv
			   > (rotate_colors ? 190.0 : 100.0) )
			    erase = TRUE;
		    }
		    else
		    {
			lost_stars = sd.live_erase - sp.live_stars;
			
			if( sd.erase_disps >
			   fv*( (rotate_colors ? 3800.0 : 2000.0)
			       - (1200*lost_stars*lost_stars)
			       / (sp.live_stars+2.0)
			       )
			   )
			    erase = TRUE;
		    }
	    }
	
	
	    /* if we can't see any stars, then we should start over */
	    if( sd.num_seen == 0 )
	    {
		if( verbose > 1 )
		    fprintf( stderr, "%s:%d  no stars are visible...\n", __FILE__, __LINE__ );
		
		new_star( prev, cur, m, vk, ak, max_stars, &sp, &sd, 0 );
	    }
	    sd.num_seen = 0;
	
	
	    /*
	     * check for X events and user commands
	     */

	    do
	    {
		if( check_event )
		{
		    XNextEvent(display.dpy, &xev);
		    HandleEvent(&xev, &sd);
		}
		
		/* Add a collapsar or a star */
		if( add || add_star )
		{
		    found = -1;
		    for( i = 0; i < max_stars; i++ )
			if( m[i] <= 0 )
			{
			    found = i;
			    break;
			}
		    
		    if( (found == -1 || stars + MAX_COLLAPSAR >= max_stars )
		       && stars + MAX_COLLAPSAR <= MAX_NUM_STARS )
		    {
			max_stars++;
			if( !num_disp_pt && sd.max_points < max_stars*10 )
			    sd.max_points = max_stars*10;
			
			
			sd.points = (XPoint *) Safe_Realloc(sd.points, sizeof(XPoint) * sd.max_points );
			sd.star_color = (int *)Safe_Realloc( sd.star_color, sizeof( int ) * max_stars );
			if( !num_disp_pt )
			{
			    sd.tmp_pts = (XPoint *) Safe_Realloc(sd.tmp_pts, sizeof(XPoint) * sd.max_points );
			    sd.pixels = (int *) Safe_Realloc(sd.pixels, sizeof(int) * sd.max_points);
			}
			cur = (point_2d *) Safe_Realloc(cur, sizeof(point_2d) * max_stars);
			prev = (point_2d *) Safe_Realloc(prev, sizeof(point_2d) * max_stars);
			sd.last_disp = (XPoint *) Safe_Realloc( sd.last_disp, sizeof(XPoint) * max_stars );
			m = (double *) Safe_Realloc(m, sizeof(double) * max_stars);
			for( k = 0; k < K; k++ )
			{
			    vk[k] = (point_2d *)Safe_Realloc( vk[k], sizeof( point_2d ) * max_stars );
			    ak[k] = (point_2d *)Safe_Realloc( ak[k], sizeof( point_2d ) * max_stars );
			}
			
			found = max_stars - 1;
			m[found] = 0;
			sd.star_color[ found ] = sd.rnd_colors[ found % NUM_COLORS ];
		    }
		    
		    if( found != -1 )
		    {
			sp.num_add++;
			stars++;
			
			new_star( prev, cur, m, vk, ak, max_stars, &sp, &sd,
				 add );
		    }
		    
		    
		    add = FALSE;
		    add_star = FALSE;
		    
		}
		
		/* delete a star */
		if( del_star )
		{
		    del_star = FALSE;
		    
		    /* kill the least massive one */
		    found_mass = collapsar*10;
		    for( i = 0, found = -1; i < max_stars; i++ )
		    {
			dist = sqrt(cur[i].x*cur[i].x + cur[i].y*cur[i].y);
			if( dist < 50 )
			    dist = 50;
			
			tmass = m[i] / dist;
			
			if( m[i] > 0 && tmass < found_mass )
			{
			    found = i;
			    found_mass = tmass;
			}
		    }
		    
		    if( found != -1 && stars > 2 )
		    {
			stars--;
			sp.live_stars--;
			sd.tstamp = time(NULL);
			set_buffer_factor( &sd, &sp );
			if( verbose > 1 )
			    fprintf( stderr, "%s:%d  live_stars=%d  m[%d]=%g\n", __FILE__, __LINE__, sp.live_stars, found, m[found] );
			
			m[found] = 0.0;
			cur[found].x = cur[found].y = prev[found].x = prev[found].y = -1.0;
			vk[0][found].x = vk[0][found].y = 0.0;
			
			if( sp.live_stars <= sp.min_stars )
			    init_stars( prev, cur, m, vk, ak, max_stars, &sp, &sd );
			
			init_fna( max_stars, m, cur, vk[0] );
		    }
		}
		
		
		/* Erase trails */
		if( erase )
		{
		    erase = FALSE;
		    
		    update_screen( &sp, &sd, &display );
		    XFlush( display.dpy );
		    sd.points_plotted += sd.points_used;
		    sd.total_points += sd.points_plotted;
		    num_erase++;
		    sd.points_plotted = 0;
		    sd.points_used = 0;
		    
		    XFillRectangle(display.dpy, display.win, display.erase_gc,
				   0,0, winW, winH);
		    
		    sd.erase_disps = 0;
		    sd.live_erase = sp.live_stars;
		    
		    clear_disp_pt( &sd );
		    plot_collapsars( cur, m, &sp, &sd, &display );
		}
		
		
		/* new start trails */
		if( new_start )
		{
		    new_start = FALSE;
		    init_stars( prev, cur, m, vk, ak, max_stars, &sp, &sd );
		}
		
		
		/* stop updating */
		if( pause_updt )
		{
		    pause_updt = FALSE;
		    paused = !paused;
		}
		
		
		/* Clean up and shut down. */
		if( stop )
		{
		    stop = FALSE; /* reset the "stop" variable */
		    
		    
		    XFillRectangle(display.dpy, display.win, display.erase_gc,
				   0,0, winW, winH);
		    XFlush( display.dpy );
		    
		    if( !root )
			XDestroyWindow(display.dpy, display.win);
		    
		    /* Free Memory */
		    free(sd.points);
		    free(sd.star_color);
		    if( !num_disp_pt )
		    {
			free(sd.tmp_pts);
			free(sd.pixels);
		    }
		    free(cur);
		    free(prev);
		    free(sd.last_disp);
		    free(m);
		    for( k = 0; k < K; k++ )
		    {
			free( vk[k] );
			free( ak[k] );
		    }
		    
		    /* Free the graphics contexts. */
		    XFreeGC(display.dpy, display.star_gc);
		    XFreeGC(display.dpy, display.erase_gc);
		    
		    if( verbose > 0 && last_init != 0 )
		    {
			sd.tstamp = time(NULL);
			if( verbose == 1 )
			    fprintf( stderr,
				    "%5d/%2d/%2d %4d %8.2f %2d %3d %3d %5d   | %4ld %5d %7d %4d%9d\n",
				    orig_sp.live_stars, orig_sp.live_stars-sp.mono_stars,
				    orig_sp.num_add, sp.num_collapsar,
				    sp.size, sp.star_circle, sp.star_line, sp.do_bounce,
				    sp.no_speed,
				    sd.tstamp - last_init, num_collide, num_edge,
				    sp.num_add, sd.num_steps );
			else
			    fprintf( stderr, "%s:%d  total life span of the star system: %ld sec (%d steps)\n", __FILE__, __LINE__, sd.tstamp - last_init, sd.num_steps );
			
			if( verbose > 3 )
			{
			    prt_sys_const( max_stars, cur, m, vk[0] );
			    prt_sys_const_delta( orig_sc, max_stars, cur, m, vk[0] );
			}
		    }
		    
		    return;
		}
		
		if( screen_exposed )
		{
		    screen_exposed = FALSE;
		    
		    redraw_screen( cur, m, &sp, &sd, &display );

		    set_redraw_full( &sd );
		}
		
		if( do_redraw )
		{
		    do_redraw = FALSE;
		    
		    XFillRectangle(display.dpy, display.win, display.erase_gc,
				   0,0, winW, winH);
		    
		    redraw_screen( cur, m, &sp, &sd, &display );
		}
		
		if( check_event )
		    check_event = XEventsQueued( display.dpy, QueuedAfterReading );
	    }
	    while( check_event );
	}
    }
}
