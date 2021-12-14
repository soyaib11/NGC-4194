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

/* note:  new_star() has similar functionality... change both places */

/*
 * The routine set_xmva sets up the star system's location (x), mass
 * (m), velocity (v) and acceleration (a).  The star system is
 * carefully crafted so that the stars orbit each other and the star
 * system doesn't drift off the screen.
 *
 *
 * The algorithm that sets up the star system is somewhat crude.  It
 * picks a direction that the star will move by combining the vector
 * needed to orbit the center and a vector needed to orbit the nearest
 * star.  The magnitude of the velocity is set by assuming that having
 * the kinetic energy to be no more than half the binding energy is a
 * good value.  This seems to work fairly well, but it still has too
 * much voodoo magic in it for my liking.
 *
 * Over the course of developing xstar, I tried many different methods
 * of setting up the initial star system.  Xgrav just set random
 * velocities and it worked well enough to be interesting.  My first
 * attempt at an improvement was to use all sorts of magic numbers and
 * weird ad-hoc formulas to try and make these random settings more
 * useful.  It worked better than xgrav's method, but it didn't scale
 * well and it clearly had room for improvement.  I found that just
 * making every star orbit the center was a major improvement, but
 * many other tweaks failed.
 *
 * I tried setting up the stars so that the velocity of each star was
 * the sum of the velocities needed to orbit each star in the star
 * system.  I had thought that would make a "perfect" orbit, but I was
 * mistaken.  The resulting velocities were too large.  The velocities
 * could be fudged to ok values, but that didn't solve all the the
 * problems.  The resulting star system tended to not orbit the
 * center, and there were still a lot of early collisions.  This was
 * before I had the binding energy calculation working, so maybe with
 * this correction factor this method might work better.
 *
 * My next though was that for each star, the remaining stars set up
 * a force field, and the gradient of that field would tell me which
 * direction the change in acceleration would be the greatest.  By picking
 * a velocity perpendicular to the gradient, the star would be moving
 * in a direction of constant acceleration and you could use the formula
 * a = v^2 / r to determine the speed.  I still don't think this would
 * create the perfect orbit for the star since the force field changes
 * with time.
 *
 * I tried to think of some sort of iterative approach to the problem,
 * where I would pick an initial approximate velocity vector for each
 * star and then let the velocities converge toward a perfect a perfect
 * orbit, but I couldn't think of any way of doing this.
 *
 * So, while the current star system is a hack, it is the best hack
 * that I know how to make.
 *
 * Heck, it is possible that if you just used random velocities like xgrav
 * had, and use the binding energy correction, that the star system
 * would work fairly well.  The binding energy correction was clearly
 * a very powerful tool to making interesting star systems.
 *
 *
 * Note:  I was just thinking...  Right now we create systems that rotate,
 * yielding systems with a large angular momentum.   Since angular momentum
 * is conserved, that means after a bunch of collisions, we must still have
 * a large angular momentum therefore the star system must be very spread
 * out.  If we want to keep things on the screen, we should try to keep the
 * angular momentum close to zero.  Or, am I missing something...
 */

double sqr( double a )
{
    return a*a;
}


void set_xmva( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp )
{
    int		i, j;

    point_2d	tp;			/* total momentum		*/
    point_2d	cm;			/* center of mass		*/
    double	tm;			/* total mass			*/

    double	dx, dy;

    double	dist;
    double	max_dist, cur_dist;

    int		cur_dist_idx;
    int		tot_stars;


    double	theta, theta0;

    point_2d	*v_0 = vk[0];

    point_2d	drift_amount;

    cm.x = cm.y = 0;
    tm = 0;
    
    
    /*
     * calculate the amount to drift
     */
    if( sp->drift )
    {
	if( sp->star_line )
	{
	    drift_amount.x = 0;
	    do
		drift_amount.y = (drand48() - .5 ) * fv_inv * .0025;
	    while( fabs( drift_amount.y ) < fv_inv * .0005 );
	}
	else
	{
	    drift_amount.x = (drand48() - .5 ) * fv_inv * .001;
	    drift_amount.y = (drand48() - .5 ) * fv_inv * .001;
	}
    }
    

    /*
     * create the star system
     */
    if( sp->star_line )
    {
	/*
	 * create stars in a straight line
	 */
	
	cm.x = cm.y = 0;
	tm = 0;

	for( i = 0; i < sp->live_stars; i++ )
	{
	    if( sp->rnd_spacing )
		prev[i].x = cur[i].x = sp->cir_dist * ( i + drand48() - .5 );
	    else
		prev[i].x = cur[i].x =
		    sp->cir_dist * (i - sp->live_stars*.5 + .5);

	    if( sp->drift )
		prev[i].y = cur[i].y = 0;
	    else
	    {
		if( i*2 < sp->live_stars )
		    prev[i].y = cur[i].y = -20;
		else
		    prev[i].y = cur[i].y = +20;
	    }


	    m[i] = fm*1;
	

	    /* keep track of what we have done for later adjustments */
	    tm += m[i];

	    cm.x += m[i]*prev[i].x;
	    cm.y += m[i]*prev[i].y;
	}
    }
	
    else if( sp->star_circle )
    {
	/*
	 * create stars in a circle
	 */

	cm.x = cm.y = 0;
	tm = 0;

	theta = 2 * M_PI/ sp->live_stars;
	theta0 = M_PI * drand48();

	for( i = 0; i < sp->live_stars; i++ )
	{
	    prev[i].x = cur[i].x = sp->cir_dist * cos( i * theta + theta0 );
	    prev[i].y = cur[i].y = sp->cir_dist * sin( i * theta + theta0 );

	    m[i] = fm*1;

	    /* keep track of what we have done for later adjustments */
	    tm += m[i];

	    cm.x += m[i]*prev[i].x;
	    cm.y += m[i]*prev[i].y;
	}
    }

    else
    {
	/*
	 * create a bunch of stars
	 */
	
	max_dist = far_dist * far_dist;
	for( i = 0; i < sp->live_stars; i++ )
	{
	    if( i < sp->mono_stars )
	    {
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
		    prev[i].x = cur[i].x = dist * cos( theta );
		    prev[i].y = cur[i].y = dist * sin( theta );

		    
		    /*
		     * don't put stars too close together
		     */
		    cur_dist = far_dist * far_dist;
		    for( j = 0; j < i; j++ )
		    {
			dx = prev[i].x - prev[j].x;
			dy = prev[i].y - prev[j].y;
			dist = dx*dx + dy*dy;
			
			if( dist < cur_dist )
			    cur_dist = dist;
		    }
		}
		while( cur_dist < 30*30
		      && ( cur_dist < 20*20 || drand48() < .5 )
		      && ( cur_dist < 10*10 || drand48() < .75 )
		      );
		
		m[i] = fm*90.+fm*RAND(120);
	    }
	    else
	    {
		/*
		 * create binary stars on the outside
		 */
		cur_dist = -1;
		cur_dist_idx = -1;
		for( j = 0; j < sp->mono_stars; j++ )
		{
		    dist = prev[j].x*prev[j].x + prev[j].y*prev[j].y;
		    if( dist < max_dist && dist > cur_dist )
		    {
			cur_dist = dist;
			cur_dist_idx = j;
		    }
		}
		max_dist = cur_dist * .9999;
		
		theta = drand48() * M_PI * 2;
		dist = 6 + drand48() * 8;
		
		prev[i].x = cur[i].x = cur[cur_dist_idx].x + dist * cos(theta);
		prev[i].y = cur[i].y = cur[cur_dist_idx].y + dist * sin(theta);
		
		
		/* cut the mass down on binary stars */
		m[i] = (fm*90.+fm*RAND(120))*.5;
		
		m[cur_dist_idx] *= .5;
		
		tm -= m[cur_dist_idx];
		cm.x -= m[cur_dist_idx]*prev[cur_dist_idx].x;
		cm.y -= m[cur_dist_idx]*prev[cur_dist_idx].y;
		
		
		/* make only half of the outer most stars binary stars */
		cur_dist = -1;
		for( j = 0; j < sp->mono_stars; j++ )
		{
		    dist = prev[j].x*prev[j].x + prev[j].y*prev[j].y;
		    if( dist < max_dist && dist > cur_dist )
			cur_dist = dist;
		}
		max_dist = cur_dist * .9999;
	    }
	    
	    
	    /* keep track of what we have done for later adjustments */
	    tm += m[i];
	    
	    cm.x += m[i]*prev[i].x;
	    cm.y += m[i]*prev[i].y;
	}
    }


    /*
     * relocate the center of mass to the origin
     */
    cm.x /= tm;
    cm.y /= tm;

    if( fabs( cm.x ) < 1 )
	cm.x = 0;
    if( fabs( cm.y ) < 1 )
	cm.y = 0;

    for( i = 0; i < sp->live_stars; i++ )
    {
	prev[i].x = cur[i].x = prev[i].x - cm.x;
	prev[i].y = cur[i].y = prev[i].y - cm.y;
    }


    /*
     * scale down the total mass by the number of stars
     *
     * I am really not sure what should be done here.
     *
     * On one hand, if you don't take into account how many stars there
     * are, and you make the average star a certain mass, then after
     * several (many?) collisions, you can end up with a very massive
     * star and you will run into a lot of overstep problems.  You
     * will also be increasing the total mass of the system as the
     * number of stars increases and that will make stars move faster
     * and thus change the accuracy of the system.
     *
     * On the other hand, if you make the star system have a constant
     * value and adjust the mass of each star accordingly, then for
     * very small star systems, the individual stars will be quite
     * massive and you will run into the overstep problem.  Also, even
     * though the total mass of the star system is constant, with lots
     * of stars the mass will be very spread out and the stars will
     * end up moving slower, thus changing the accuracy of the system.
     *
     * So, the mass of the star system should scale with the number of
     * stars, but sub-linearly.  what this formula is, I don't know
     * and I don't know how to figure it out.
     *
     * For right now, it appears that keeping a constant star system
     * mass is better than a simple linear scale, so that's what I am
     * doing.
     */
    if( sp->num_collapsar || sp->do_bounce )
	tm = (fm * 4) / tm;
    else if( sp->few_stars && sp->live_stars > 2 )
	tm = (fm * 6) / tm;
    else
	tm = (fm * 16) / tm;
    
	
    if( tm <= 0 )
    {
	fprintf( stderr, "%s:%d  Error!  tm=%g  should be >0\n",
		__FILE__, __LINE__, tm );
	tm = fm*.1;
    }

    for( i = 0; i < sp->live_stars; i++ )
    {
	m[i] *= tm;
	if( m[i] >= collapsar )
	{
	    /* This _should_, by construction, never happen. */
	    if( verbose > 0 )
		fprintf( stderr, "%s:%d  scaled mass too large: m[%d]=%g\n",
			__FILE__, __LINE__, i, m[i]/fm );

	    m[i] = collapsar * .999;
	}
    }


    /*
     * create the collapsars
     */
    tot_stars = sp->live_stars + sp->num_collapsar;
    if( sp->num_collapsar )
    {
	cm.x = cm.y = 0;
	tm = 0;
	
	for( i = sp->live_stars; i < tot_stars; i++ )
	{
	    /*
	     * don't put stars too close together and use the appropriate
	     * distribution of stars
	     */
	    do
	    {
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
		    dist = 100;
		}
		if( !sp->do_bounce )
		{
		    prev[i].x = cur[i].x = dist * .5 * cos( theta );
		    prev[i].y = cur[i].y = dist * .5 * sin( theta );
		}
		else
		{
		    prev[i].x = cur[i].x = dist * cos( theta );
		    prev[i].y = cur[i].y = dist * sin( theta );
		}

		cur_dist = far_dist * far_dist;
		for( j = sp->live_stars; j < i; j++ )
		{
		    dx = prev[i].x - prev[j].x;
		    dy = prev[i].y - prev[j].y;
		    dist = dx*dx + dy*dy;

		    if( dist < cur_dist )
			cur_dist = dist;
		}
	    }
	    while( cur_dist < 30*30
		  && ( cur_dist < 20*20 || drand48() < .5 )
		  && ( cur_dist < 10*10 || drand48() < .75 )
		  );

	    v_0[i].x = 0.0;
	    v_0[i].y = 0.0;
	    m[i] = collapsar * 1.0001;
	    
	    /* keep track of what we have done for later adjustments */
	    tm += m[i];
	    
	    cm.x += m[i]*prev[i].x;
	    cm.y += m[i]*prev[i].y;
	}
	
	
	/*
	 * relocate the collapsars' center of mass to the origin
	 */
	cm.x /= tm;
	cm.y /= tm;
	
	for( i = sp->live_stars; i < tot_stars; i++ )
	{
	    prev[i].x = cur[i].x = prev[i].x - cm.x;
	    prev[i].y = cur[i].y = prev[i].y - cm.y;
	}
    }


    /*
     * set the velocities
     */
    for( i = 0; i < sp->live_stars; i++ )
	set_va( cur, prev, m, vk, ak, sp, i );
    


    /*
     * make sure the system has no overall momentum
     */
    tm = 0;
    tp.x = tp.y = 0;
    for( i = 0; i < sp->live_stars; i++ )
    {
	tp.x += v_0[i].x*m[i];
	tp.y += v_0[i].y*m[i];

	tm += m[i];
    }

    tp.x /= tm;
    tp.y /= tm;

    /* ok, if there aren't many stars, let them drift around */
    if( sp->drift )
    {
	tp.x += drift_amount.x;
	tp.y += drift_amount.y;
    }

    for( i = 0; i < sp->live_stars; i++ )
    {
	v_0[i].x -= tp.x;
	v_0[i].y -= tp.y;
    }


    /*
     * if we are drifting, move the starting locations of the stars
     * in the opposite direction
     */

    if( sp->drift && ( !sp->do_bounce || sp->star_line ) )
    {
	double		gap;
	point_2d	edge;
	point_2d	offset;

	if( sp->star_line )
	    gap = 40;
	else
	    gap = sp->size * 1.5;
	    
	edge.x = ( max_x < -min_x ) ? max_x : -min_x;
	edge.y = ( max_y < -min_y ) ? max_y : -min_y;

	if( gap > edge.x )
	    gap = edge.x;
	if( gap > edge.y )
	    gap = edge.y;

	if( fabs( drift_amount.x / edge.x ) > fabs( drift_amount.y / edge.y ) )
	{
	    /* the x direction is the limiting factor */
	    offset.x = ( drift_amount.x > 0 ) ? -edge.x + gap : edge.x - gap;
	    offset.y = offset.x * drift_amount.y / drift_amount.x;
	}
	else
	{
	    /* the y direction is the limiting factor */
	    offset.y = ( drift_amount.y > 0 ) ? -edge.y + gap : edge.y - gap;
	    offset.x = offset.y * drift_amount.x / drift_amount.y;
	}

	for( i = 0; i < tot_stars; i++ )
	{
	    prev[i].x = cur[i].x = prev[i].x - offset.x;
	    prev[i].y = cur[i].y = prev[i].y - offset.y;
	}

    }
}


void set_va( point_2d *cur, point_2d *prev, double *m,
	    point_2d *vk[], point_2d *ak[],
	    sys_param_type *sp,
	    int i )
{
    int		j;

    double	tm;			/* total mass			*/

    double	dx, dy;
    double	tdx, tdy;

    double	dist, tdist;

    double	k, u;
    double	v_fac;

    int		tot_stars;


    point_2d	*v_0 = vk[0];
    point_2d	*a_0 = ak[0];

    double	max_energy;

    double	a_mag, q_a_mag;
    double	f;


    tot_stars = sp->live_stars + sp->num_collapsar;

    
    /*
     * initialize the acceleration array
     */
    for( j = 0; j < tot_stars; j++ )
    {
	double f, f_d;

	if( j == i )
	    continue;

	if( m[j] <= 0.0 )
	    continue;
	    
	dx = prev[j].x - prev[i].x;
	dy = prev[j].y - prev[i].y;
	dist = dx*dx + dy*dy;
	    
	if( dist >= collide_dist )
	{
	    f = G / (sqrt(dist) * dist);
		
	    f_d = dx * f;		
	    a_0[i].x += f_d * m[j];
		
	    f_d = dy * f;
	    a_0[i].y += f_d * m[j];
	}
    }
    
    
    /*
     * set the velocities
     */
    dx = cur[i].x;
    dy = cur[i].y;
    dist = sqrt( dx*dx + dy*dy );
    
    if( dist > 4 && !sp->star_circle )
    {
	/*
	 * assume that the stars form a uniformly dense sphere and
	 * calculate the gravitational force on this particular star.
	 *
	 * we assume that all stars closer to the center than this one
	 * are uniformly distributed and thus can be treated as a
	 * point mass at the center.  we also assume that the stars
	 * further out than this one are uniformly distributed and thus
	 * have no affect on this star.  Of course, none of this is true.
	 */
	tm = 0;
	dist *= .99999999;
	for( j = 0; j < tot_stars; j++ )
	{
	    if( m[j] <= 0 )
		continue;
	    
	    tdx = cur[j].x;
	    tdy = cur[j].y;
	    tdist = sqrt( tdx*tdx + tdy*tdy );
	    
	    if( tdist < dist )
		tm += m[j];
	}
	if( tm < .5*fm )
	    tm = .5*fm;
	
	if( sp->num_collapsar )
	    tm *= .30;
	else
	    tm *= .45;
	
	
	/* vector tangent to a perfectly circular orbit */
	f = sqrt( G * tm / dist);
	
	v_0[i].x = (dy / dist) * f;
	v_0[i].y = (-dx / dist) * f;
	
    }
    
    
    /*
     * add a twist to the velocity to avoid any close by stars
     */
    a_mag = sqrt( a_0[i].x * a_0[i].x + a_0[i].y * a_0[i].y );
    if( a_mag*fv*fv > 1E-9 )
    {
	q_a_mag = sqrt( a_mag );		
	
	if( sp->star_circle )
	{
	    /* orbit center of the stars system */
	    v_0[i].x += (-a_0[i].y * sqrt(sp->cir_dist)) / q_a_mag;
	    v_0[i].y += (a_0[i].x * sqrt(sp->cir_dist)) / q_a_mag;
	}
	else
	{
	    dist = far_dist*far_dist;
	    for( j = 0; j < tot_stars; j++ )
	    {
		if( i == j )
		    continue;

		if( m[j] <= 0 )
		    continue;
		
		tdx = (prev[i].x*m[i] + prev[j].x*m[j]) / (m[i] + m[j])
		    - prev[i].x;
		tdy = (prev[i].y*m[i] + prev[j].y*m[j]) / (m[i] + m[j])
		    - prev[i].y;
		
		tdist = tdx*tdx + tdy*tdy;
		
		if( tdist < dist )
		    dist = tdist;
	    }
	    dist = sqrt( sqrt( dist ) );


	    /* orbit closest stars with radius of 'dist' to the C of M */
	    /* derived from a = v^2/rcm and  a = r12/|r12| * G * m1 / r12^2 */
	    if( sp->min_angular_mom )
	    {
		v_0[i].x += (a_0[i].y * 6 * dist) / q_a_mag;
		v_0[i].y += (-a_0[i].x * 6 * dist) / q_a_mag;
	    }
	    else
	    {
		v_0[i].x += (-a_0[i].y * dist) / q_a_mag;
		v_0[i].y += (a_0[i].x * dist) / q_a_mag;
	    }


	    /* give systems with collapsars a random twist */
	    if( sp->num_collapsar )
	    {
		v_0[i].x *= 1 + (drand48() - .4)/2;
		v_0[i].y *= 1 + (drand48() - .4)/2;
	    }
	}
    }
    
    
    /*
     * make sure this star is still bound
     */
    if( !sp->no_speed )
    {
	k = kinetic_energy( i, tot_stars, cur, m, v_0 );
	u = binding_energy( i, tot_stars, cur, m, v_0 );
	
	if( sp->calc_energy_fact )
	{
	    if( sp->num_collapsar )
	    {
		sp->energy_fact = drand48() * .25 + .30;
	    }
	    else
		switch( sp->star_distrib )
		{
		case DIST_CENTER:
		    sp->energy_fact = drand48() * .1 + .50;
		    break;
		    
		case DIST_UNIFORM:
		    sp->energy_fact = drand48() * .1 + .45;
		    break;
		    
		case DIST_RING:
		    sp->energy_fact = drand48() * .1 + .40;
		    break;
		    
		default:
		    sp->energy_fact = 0;
		    fprintf( stderr, "%s:%d  Error!  star distribution=%d\n",
			    __FILE__, __LINE__, sp->star_distrib );
		    exit( 1 );
		    break;
		}
	}
	max_energy = -u * sp->energy_fact;
	
	
	if( k > max_energy || sp->star_circle )
	{
	    v_fac = sqrt( max_energy / k );
	    
	    v_0[i].x *= v_fac;
	    v_0[i].y *= v_fac;
	}
#if 0
	else
	    v_fac = -1;
	    
	    k = kinetic_energy( i, tot_stars, cur, m, v_0 );
	    u = binding_energy( i, tot_stars, cur, m, v_0 );
	    
	    fprintf( stderr, "%s:%d %d  k=%9.6f u=%9.6f  -k/u=%9.6f   |%7.4f,%7.4f|  v_fac=%g\n", __FILE__, __LINE__, i, k*fv*fv/fm, u*fv*fv/fm, -k/u, fv*v_0[i].x, fv*v_0[i].y, v_fac );
#endif
	    
    }
    
    
    /* sometimes it's nice to star with no velocity */
    if( sp->no_speed )
	v_0[i].x = v_0[i].y = 0;
}
