#include "xstar.h"
#include "xstar_ext.h"


/*
 * Check_bounce() is _really_ expensive.  Many banks charge $25 for the first
 * time, and many states require you to take 'classes' if you bounce checks.
 */

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
 * This code checks to see if any stars should bounce off of each other.  It
 * mostly works, but for multi-step routines like taylor3 or the Adam-Bashford
 * routines, it is fairly slow because it requires a movement restart.  There
 * are other inefficiencies in the code as well.  Fortunately, it only
 * calculates these expensive formulas when an actual collision is happening,
 * and that doesn't happen very often.
 *
 * This routine does not work in all cases though.  For some reason, if you
 * have more than two stars colliding at the same spot, the results are not
 * correct.  You can get this when you take a star circle, enable the bounce
 * code and start with no speed.  Then, all the stars drop toward the center
 * of the circle and when they get there, they should all bounce back to the
 * edge along the same path.  This only happens for n=2,4, for other values
 * n, the stars go off at random angles.
 *
 * [ I have checked into this a little bit more... It appears that the bounce
 *   code is working "correctly", if you take the "duration of the impulse" as
 *   being zero (or close too it) and you count stars with lower index numbers
 *   to be "slightly" ahead of the stars with higher index numbers.  Besides
 *   the fix described in the next paragraph, this problem might also be
 *   "fixed" by accumulating all the bounce effects and then somehow applying
 *   them after all stars have been checked.  I am not sure that this would
 *   be simpler than the following, better, fix... ]
 *
 * I have thought about how maybe we should bounce stars via the same method
 * that real things bounce.  That is, when objects get close enough together,
 * the atoms no longer look like neutral balls, and instead the electrical
 * repulsion between the electrons stars to take over.  If we created a very
 * short distance, very large magnitude negative force field around each
 * star, I think we could get rid of the above problem.  This would create
 * the non-zero duration of the impulse that real world items have.
 */

#define bounce_dist	(5)		/* large enough to be visible on scr */

/* #define	bounce_debug */


void check_bounce( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd )
{
    register double	m_i, m_j;
    
    register int        i, j;		/* star index */
    
    point_2d	*v_1 = vk[1];
    
    
    
    /* bounce stars off of each other */
    for( i = 0; i < max_stars; i++ )
    {
	m_i = m[i];
	if( m_i <= 0.0 )
	    continue;
	
	
	for( j = i+1; j < max_stars; j++ )
	{
	    double dx, dy, n2;
	    point_2d tvi, tvj;
	    
	    m_j = m[j];
	    if( m_j <= 0.0 )
		continue;
	    
	    dx = prev[j].x - prev[i].x;
	    dy = prev[j].y - prev[i].y;
	    n2 = dx*dx + dy*dy;
	    
	    
	    /* this is probably too large of a value... */
	    if( (m[j] < collapsar && n2 < bounce_dist*bounce_dist )
	       || (m[j] >= collapsar
		   && n2 < bounce_dist*bounce_dist*bounce_dist
		   )
	       )
	    {
		/* collision:  perfectly elastic -> bounce */
		
		double n = sqrt(n2);
		
		
		/* projected velocities along alternate coordinate sys */
		point_2d	v1xp, v1yp, v2xp, v2yp;
		point_2d	xp, yp;
		
		/*
		 * figure out new coordinate system
		 *
		 * The x-prime axis is along the line connecting the
		 * two bodies.  The y-prime axis is "flat wall" that the
		 * bodies appear to bounce off of.  The momentum and
		 * velocity along the y-prime axis is not changed by
		 * the collision.
		 */
		xp.x = dx / n;
		xp.y = dy / n;
		
		yp.x = xp.y;
		yp.y = -xp.x;

		
		/* figure out projections of velocities along new axes */
		v1yp.x = (yp.x * v_1[i].x + yp.y * v_1[i].y) * yp.x;
		v1yp.y = (yp.x * v_1[i].x + yp.y * v_1[i].y) * yp.y;
		
		v2yp.x = (yp.x * v_1[j].x + yp.y * v_1[j].y) * yp.x;
		v2yp.y = (yp.x * v_1[j].x + yp.y * v_1[j].y) * yp.y;
		
		v1xp.x = (xp.x * v_1[i].x + xp.y * v_1[i].y) * xp.x;
		v1xp.y = (xp.x * v_1[i].x + xp.y * v_1[i].y) * xp.y;
		
		v2xp.x = (xp.x * v_1[j].x + xp.y * v_1[j].y) * xp.x;
		v2xp.y = (xp.x * v_1[j].x + xp.y * v_1[j].y) * xp.y;
		
		/*
		 * basically, we just need to treat the collision as a
		 * one dimensional collision along the x-prime axis.
		 */
		if( m_i >= collapsar )
		{
		    tvi = v_1[i];
		
		    tvj.x = v2yp.x - v2xp.x;
		    tvj.y = v2yp.y - v2xp.y;
		}
		else if( m_j >= collapsar )
		{
		    tvi.x = v1yp.x - v1xp.x;
		    tvi.y = v1yp.y - v1xp.y;

		    tvj = v_1[j];
		}
		else
		{
		    double  tm_inv = 1. / (m_i + m_j);
		    tvi.x = v1yp.x + (m_i * v1xp.x - m_j * v1xp.x + 2 * m_j * v2xp.x) * tm_inv;
		    tvi.y = v1yp.y + (m_i * v1xp.y - m_j * v1xp.y + 2 * m_j * v2xp.y) * tm_inv;
		
		    tvj.x = v2yp.x + (m_j * v2xp.x - m_i * v2xp.x + 2 * m_i * v1xp.x) * tm_inv;
		    tvj.y = v2yp.y + (m_j * v2xp.y - m_i * v2xp.y + 2 * m_i * v1xp.y) * tm_inv;
		}
		
		
		/*
		 * There has got to be a better way to see if the
		 * stars are moving away from each other than
		 * checking the distances like this.  This is way
		 * too expensive to calculate.
		 *
		 * the problem that I am fighting is that once inside the
		 * bounce distance, this code always says that the stars
		 * are colliding, even though after the first bounce, they
		 * are moving away from each-other.
		 *
		 * we know that the tvi and tvj vectors are collinear,
		 * and that should make things easier to calculate...
		 */
		dx = (prev[j].x + tvj.x) - (prev[i].x + tvi.x);
		dy = (prev[j].y + tvj.y) - (prev[i].y + tvi.y);
		
		if( n2 < (dx*dx + dy*dy) )
		{
#ifdef bounce_debug
		    static double	last_mag = 0, cur_mag = 0;
		    
		    fprintf( stderr, "\n%s:%d  i=(%g,%g)  j=(%g,%g)  delta=(%g,%g)  n=%g\n", __FILE__, __LINE__, prev[i].x, prev[i].y, prev[j].x, prev[j].y, dx, dy, n );
#define rtd(x)	((x)*360/(M_PI*2))
		    fprintf( stderr, "%s:%d  xp=(%g,%g) |%g| <%g\n",
			    __FILE__, __LINE__,
			    xp.x, xp.y,
			    sqrt((xp.x*xp.x) + (xp.y*xp.y)),
			    rtd( atan2( xp.y, xp.x ) ) );

		    fprintf( stderr, "%s:%d  yp=(%g,%g) |%g| <%g\n",
			    __FILE__, __LINE__,
			    yp.x, yp.y,
			    sqrt((yp.x*yp.x) + (yp.y*yp.y)),
			    rtd( atan2( yp.y, yp.x ) ) );

		    fprintf( stderr, "%s:%d  stars: %d,%d  v1=(%g,%g) |%g| <%g  v1yp=(%g,%g) v1xp=(%g,%g)  (%g,%g)\n",
			    __FILE__, __LINE__, i, j,
			    v_1[i].x, v_1[i].y,
			    sqrt((v_1[i].x*v_1[i].x) + (v_1[i].y*v_1[i].y)),
			    rtd( atan2( -v_1[i].y, -v_1[i].x ) ),
			    v1yp.x, v1yp.y, v1xp.x, v1xp.y,
			    v1xp.x + v1yp.x, v1yp.y + v1xp.y );
		    fprintf( stderr, "%s:%d  v1=(%g,%g) |%g| <%g\n",
			    __FILE__, __LINE__,
			    tvi.x, tvi.y,
			    sqrt((tvi.x*tvi.x) + (tvi.y*tvi.y)),
			    rtd( atan2( tvi.y, tvi.x ) ) );
		    fprintf( stderr, "%s:%d  stars: %d,%d  v2=(%g,%g) |%g| <%g  v2yp=(%g,%g) v2xp=(%g,%g)  (%g,%g)\n",
			    __FILE__, __LINE__, i, j,
			    v_1[j].x, v_1[j].y,
			    sqrt((v_1[j].x*v_1[j].x) + (v_1[j].y*v_1[j].y)),
			    rtd( atan2( -v_1[j].y, -v_1[j].x ) ),
			    v2yp.x, v2yp.y, v2xp.x, v2xp.y,
			    v2xp.x + v2yp.x, v2yp.y + v2xp.y );
		    fprintf( stderr, "%s:%d  v2=(%g,%g) |%g| <%g\n",
			    __FILE__, __LINE__,
			    tvj.x, tvj.y,
			    sqrt((tvj.x*tvj.x) + (tvj.y*tvj.y)),
			    rtd( atan2( tvj.y, tvj.x ) ) );
#endif
		    
		    
		    v_1[i] = tvi;
		    v_1[j] = tvj;

		    restart_fna( max_stars, m, prev, v_1 );

#ifdef bounce_debug
#if 0
		    cur_mag = sqrt((v_1[i].x*v_1[i].x) + (v_1[i].y*v_1[i].y));
		    
		    fprintf( stderr, "%s:%d  v1=(%g,%g) %g  diff=%.3g%%\n",
			    __FILE__, __LINE__,
			    v_1[i].x, v_1[i].y, cur_mag,
			    100.*(last_mag/cur_mag-1.) );
		    last_mag = cur_mag;
#endif
#endif
		}
#ifdef bounce_debug
		else
		    fprintf( stderr, "%s:%d  ignoring the bounce...\n", __FILE__, __LINE__ );
#endif
	    }
	}
	
    }
}
