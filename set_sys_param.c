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
 * This routine decides what kind of star system to create.  It
 * selects from either a normal star system, or one of many special
 * kinds of star systems.  It then sets all the appropriate parameters
 * and lets set_xmva() do the real work.  (Just deciding what do
 * do seems like a fairly good sized routine...)
 */


/* #define DEBUG */


void set_sys_param( sys_param_type *sp, int max_stars )
{
    int		tstars;
#if defined( DEBUG )
    static int  seed = 1;
    long	junk;
    int		i;
#endif
    
#if defined( DEBUG )
    /*
     * set up the random number generator to a predictable state...
     */
    srand48( seed++ );

    /* seeds won't change much, so use the rng a few times to let the
     * sequences diverge a little bit... */
    for( i = 0; i < 48; i++ )
	junk = lrand48();
#endif
    
#if 0
/* ** debug ** */
#if 0
    /* for easier debugging */
    sp->star_distrib = DIST_UNIFORM;
    sp->num_collapsar = 0;
    sp->star_circle = 1;
    sp->star_line = 0;
    sp->do_bounce = 0;
    sp->no_speed = 0;
    sp->few_stars = 0;
    sp->live_stars = stars;
    sp->cir_dist = 70;
    sp->min_stars = 1;
    sp->num_add = 0;
    sp->drift = 0;
#else
    	    /* 40% of the time, create a normal star system */
	    sp->star_distrib = lrand48() % 3;
	    sp->num_collapsar = 0;
	    sp->star_circle = 0;
	    sp->star_line = 0;
	    sp->do_bounce = 0;
	    sp->no_speed = 0;
	    sp->cir_dist = 0;
	    sp->few_stars = 0;
	    sp->live_stars = stars;
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->drift = -1;
	    sp->min_angular_mom = 1 /* lrand48() % 2 */;

#endif
#else

    sp->live_stars = 2*stars;
    sp->num_collapsar = 0;
    sp->star_circle = 0;
    sp->star_line = 0;
    sp->do_bounce = 0;
    sp->no_speed = 0;
    sp->cir_dist = 0;
    sp->few_stars = 0;
    sp->drift = -1;
    sp->min_angular_mom = 0;

    /* try creating a star system, but don't violate these restrictions */
    while( sp->live_stars > stars
	  || sp->num_collapsar + sp->live_stars > max_stars
	  || (sp->few_stars && sp->live_stars*2 > max_stars)
	  )
    {
	switch( lrand48() % 10 )
	{
	case 0:
	case 1:
	case 2:
	case 3:
	    /* 40% of the time, create a normal star system */
	    sp->star_distrib = lrand48() % 3;
	    sp->num_collapsar = 0;
	    sp->star_circle = 0;
	    sp->star_line = 0;
	    sp->do_bounce = 0;
	    sp->no_speed = 0;
	    sp->cir_dist = 0;
	    sp->few_stars = 0;
	    sp->live_stars = stars;
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->drift = -1;
	    sp->min_angular_mom = lrand48() % 2;
	    break;
		

	case 4:
	    /* 10% of the time, create a small number of stars */
	    sp->star_distrib = lrand48() % 3;
	    sp->num_collapsar = 0;
	    sp->star_circle = 0;
	    sp->star_line = 0;
	    sp->do_bounce = 0;
	    sp->no_speed = 0;
	    sp->cir_dist = 0;
	    sp->few_stars = 1;
	    sp->live_stars = 2 + lrand48() % 4;
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->drift = -1;
	    sp->min_angular_mom = 0;
	    break;
		
	case 5:
	    /* 10% of the time, create a circular star system */
	    tstars = ( stars > 28 ) ? 25 : stars-3;
	    sp->star_distrib = DIST_UNIFORM;
	    sp->num_collapsar = 0;
	    sp->star_circle = 1;
	    sp->star_line = 0;
	    sp->do_bounce = 0;
	    sp->no_speed = 0;
	    sp->few_stars = 0;
	    sp->live_stars = 3.5 + drand48() * drand48() * tstars;
	    sp->cir_dist = 50 + RAND(80);
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->drift = (lrand48() % 5) == 0;
	    sp->min_angular_mom = 0;
	    break;
		
	case 6:
	    /* 10% of the time, create a linear star system */
	    sp->star_distrib = DIST_UNIFORM;
	    sp->num_collapsar = 0;
	    sp->star_circle = 0;
	    sp->star_line = 1;
	    sp->do_bounce = 1;
	    sp->no_speed = 1;
	    sp->few_stars = 0;
	    switch( lrand48() % 8 )
	    {
	    case 0:
		sp->live_stars = 2;
		sp->drift = 1;
		break;

	    case 1:
	    case 2:
	    case 3:
		sp->live_stars = 4;
		sp->drift = lrand48() % 2;
		break;

	    default:
		tstars = ( stars > 27 ) ? 25 : stars-2;
		sp->live_stars = 5 + drand48()*drand48()*drand48() * tstars;
		sp->drift = lrand48() % 2;
		break;
	    }
	    sp->cir_dist = 40 + drand48() * 60;
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->min_angular_mom = 0;

	    /* due to bugs in the bounce code, this should be even */
	    sp->live_stars += sp->live_stars & 1;

	    /* if we only have 2 stars, put a spin on it */
	    if( sp->live_stars == 2 )
	    {
		if( (lrand48() % 4) )
		{
		    sp->no_speed = 0;
		    sp->drift = (lrand48() % 3) > 0;
		}
	    }
	    else if( sp->drift )
		sp->rnd_spacing = (lrand48() % 3) == 0;

	    break;
		
	case 7:
	case 8:
	    /* 20% of the time, create bouncing stars */
	    sp->star_distrib = lrand48() % 3;
	    sp->num_collapsar = 0;
	    sp->star_circle = 0;
	    sp->star_line = 0;
	    sp->do_bounce = 1;
	    sp->no_speed = (lrand48() % 3) == 0;
	    sp->cir_dist = 0;
	    sp->few_stars = 0;
	    sp->live_stars = 3 + lrand48() % 6;
	    sp->min_stars = 1;
	    sp->num_add = 0;
	    sp->rnd_spacing = 0;
	    sp->drift = 1;
	    sp->min_angular_mom = 0;

	    if( sp->no_speed && sp->star_distrib == DIST_CENTER )
		sp->star_distrib = DIST_UNIFORM;

	    /* every once and a while, throw in some collapsars */
	    if( lrand48() % 3 == 0 )
		sp->num_collapsar = lrand48() % (MAX_COLLAPSAR-1);
	    break;


	default:
	    /* 20% of the time, create collapsars */
		
	    sp->star_distrib = 1 + lrand48() % 2;

	    switch( lrand48() % 10 )
	    {
	    case 0:
		sp->num_collapsar = 1;
		sp->min_stars = 1;
		break;
		    
	    case 1:
	    case 2:
		sp->num_collapsar = 2;
		sp->min_stars = 0;
		break;
		    
	    case 3:
	    case 4:
	    case 5:
		sp->num_collapsar = 3;
		sp->min_stars = 0;
		break;
		    
	    case 6:
	    case 7:
		sp->num_collapsar = 4;
		sp->min_stars = 0;
		break;

	    case 8:
		sp->num_collapsar = 5;
		sp->min_stars = 0;
		break;
		    
	    case 9:
	    default:
		sp->num_collapsar = 6;
		sp->min_stars = 0;
		break;
	    }
		
	    sp->star_circle = 0;
	    sp->star_line = 0;
	    sp->do_bounce = 0;
	    sp->no_speed = 0;
	    sp->cir_dist = 0;
	    sp->few_stars = 0;
	    sp->live_stars = stars;
		
		
	    sp->num_add = stars / 3;
	    sp->live_stars -= sp->num_add;
	    sp->rnd_spacing = 0;
	    sp->drift = 0;
	    sp->min_angular_mom = 0;
	}
    }	
#endif

    /* 15-25% binary stars */
    if( sp->live_stars > 9 && !sp->star_circle && !sp->star_line )
	sp->mono_stars = sp->live_stars * (.925 - .05 * drand48())
	    + drand48()*2;
    else
	sp->mono_stars = sp->live_stars;
    if( sp->mono_stars > sp->live_stars )
	sp->mono_stars = sp->live_stars;


    /* should stars drift? */
    if( sp->drift == -1 )
	if( sp->live_stars <= 5 )
	{
	    if( sp->live_stars == 2 )
		sp->drift = (lrand48() % 3) < 2;
	    else
		sp->drift = 1;
	}
	else
	    sp->drift = 0;


    /*
     * calculate the size and amount of energy the system should have
     */
    sp->size = ( (sp->live_stars + sp->num_collapsar) * STARAREA
		+ 15000) * sqrt( sqrt( collide_dist ) );
    switch( sp->star_distrib )
    {
    case DIST_CENTER:
	sp->size = sqrt( sp->size ) * .6666;
	break;

    case DIST_UNIFORM:
	sp->size = sqrt( sp->size ) * .64;
	break;

    case DIST_RING:
	sp->size = sqrt( sp->size ) * .58;
	break;

    default:
	sp->size = 0;
	fprintf( stderr, "%s:%d  Error!  star distribution=%d\n",
		__FILE__, __LINE__, sp->star_distrib );
	exit( 1 );
	break;
    }

    if( sp->do_bounce )
	sp->size *= .5;


    /*
     * calculate the energy factor for some types of systems
     */
    sp->energy_fact = -1;
    if( sp->star_circle )
    {
	if( lrand48() % 16 )
	    sp->energy_fact = drand48() * .31 + .54;
	else
	    sp->energy_fact = .5;	/* perfect circle		*/
    }

    else if( sp->live_stars <= 2 )
    {
	if( sp->star_line )
	    sp->energy_fact = drand48() * .04 + .0005;
	else if( (lrand48() % 8) == 0 )
	    sp->energy_fact = .5;
	else if( lrand48() % 2 )
	    sp->energy_fact = drand48() * .9 + .06;
	else
	    sp->energy_fact = drand48() * .4 + .06;
    }
    sp->calc_energy_fact = ( sp->energy_fact == -1 );
    
	
    
    if( verbose > 1 )
	fprintf( stderr, "%s:%d  num stars=%d  num collapsars=%d  size=%f\n",
		__FILE__, __LINE__, sp->live_stars, sp->num_collapsar,
		sp->size );

}
