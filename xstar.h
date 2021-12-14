#ifndef include_xstar
#define include_xstar


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com) and others
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/



#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>


/* standard POSIX typedefs */
#ifdef NEEDS_TYPEDEFS
typedef	unsigned char	uchar_t;
typedef	unsigned short	ushort_t;
typedef	unsigned int	uint_t;
typedef	unsigned long	ulong_t;
#endif



/* code profiling options */

#if 0
#define MARK
#include <prof.h>
#else
#define MARK(K)	{}
#endif


/* define if you want star movement data dumped to stdout */
/* define  DUMP_STAR_DATA */


/* Constants */

/* Fudge factors
 *
 * The higher FV is, the slower the system runs, but the more accurately it
 * runs.  If you see stars shooting off very quickly when they get close
 * together, then that is a sign that FV is way too small.  If you notice that
 * the star system gains energy over time (stars are moving further out),
 * then that is a sign that FV is not large enough.  In general, the larger
 * FV, the longer the star system will continue to have "interesting" 
 * interactions because all of the stars will stay on the screen.
 *
 * If you see the user falling a sleep because things are running so
 * slowly, then that is a sign that FV is too large.  If you see both
 * things happening, that that is a sign that you need a faster
 * computer.
 *
 *
 * The higher FM is, the more massive the star systems will be and the more
 * the stars will tend to fall in toward the center.  This constant used
 * to have a major effect on init_star(), but it has since been re-written
 * and there are only a few spots where change FM will not have the same
 * effect as changing FV.
 *
 * Note:  These fudge factors are best changed by the -a command line option.
 */
#define FV	(.8)			/* larger values => slower movement */
#define FG	(16.0)			/* don't change			*/
#define FM	(1.0/(16*FG*FV*FV))	/* larger values => move attraction */

/* #define G ((FG*FG)/256.) */
#define G	(1.0)			/* gravitation const (in funny units)*/


/*
 * Units:	xstar	normal
 * distance	1	1 millimeter
 * mass		1	1 gram
 * time		1	8.16859841099e-9 s	(?? is this right?)
 *
 * G = 6.6726E-11 N*m^2/kg^2 = 6.6726E-11 (kg*m/s^2)*m^2/kg^2
 *   = 6.6726E-11 m^3/(s^2 kg)
 */

#ifndef NUM_COLORS
#define NUM_COLORS 48	  /* numbers of colors to rotate through	*/
			  /* limited to a max of 64 by disp_star_type   */
#endif


#define TRUE 1
#define FALSE 0
#define STD_STR 100
#define STARS   15			/* number of stars */
#define DELAY   0			/* microsecond delay between updates */
#define STARAREA 10000.0		/* how crowded stars should be */
#define COLLAPSAR (16.0)		/* mass of a collapsar		*/
#define WINWIDTH 512			/* default window width */
#define WINHEIGHT 512			/* default window height */
#define WINMIN	128			/* minimum window size	*/
#define ALIVE_MASK      (SubstructureNotifyMask | KeyPressMask | PointerMotionMask)
#define MAX_CVAL	65535		/* maximum value of a color	*/
#define POLL_SKIP	64		/* keep a power of 2		*/
#define MAX_COLLAPSAR	6		/* max number of initial collapsars */
#define DEFAULT_COLLIDE 1.		/* default collide distance	*/

#define MAX_NUM_STARS	((1<<8)-1)	/* maximum number of stars allowed */
					/* as limited by disp_star_type.*/
					/* the last star is used as a flag */

#define K	9			/* # of prev values of a,v to keep  */



/* Error Codes */
#define FATAL   -1
#define WARNING -2

/* Macros */
#ifdef USE_RANDOM
#define srand48(a)	srandom(a)
#define lrand48()	random()
#define drand48()	(random()/2147483647.)
#endif
#define RAND(v) ((lrand48() % v) - (v/2)) /* random number around 0 */

#define array_elem(x)	(sizeof(x)/sizeof(x[0]))
#define collide( i, j, cur, prev, m, vk, ak, sp, sd ) do_collide( i, j, cur, prev, m, vk, ak, sp, sd, __FILE__, __LINE__ )
#define prt_sys_const( n, cur, m, v_0 ) do_prt_sys_const( n, cur, m, v_0, __FILE__, __LINE__ )
#define prt_sys_const_delta( sc0, n, cur, m, v_0 ) do_prt_sys_const_delta( sc0, n, cur, m, v_0, __FILE__, __LINE__ )


/* Type Definitions */
typedef struct _disp
{
        Window  win;
        Display *dpy;
        int     screen;
        Window  root;
        char    *dname;
        long    star, bg; /* colors */
        XColor  bg_xcolor;
        GC      star_gc;
        GC      erase_gc;
        Atom    kill_atom, protocol_atom;
        Colormap cmap;
} disp;

typedef struct
{
    double	x, y;
} point_2d;


/* system constants of motion */
typedef struct
{
    double	k, u;			/* kinetic and potential energies */
    double	e;			/* total energy of the system	*/
    point_2d	cm;			/* center of mass		*/
    point_2d	p;			/* total momentum of the system	*/
    double	l;			/* tot angular momentum of sys	*/
    double	m;			/* total mass of the system	*/
} sys_const;

typedef int	dist_type;
#define	DIST_CENTER	0		/* higher probability toward center */
#define	DIST_UNIFORM	1		/* uniform probability		*/
#define	DIST_RING	2		/* higher probability toward outside */



/* system initialization parameters */
typedef struct
{
    dist_type	star_distrib;		/* distribution method for stars*/
    double	size;			/* area to spread stars over	*/
    int		num_collapsar;		/* number of collapsars to create*/
    int		mono_stars;		/* number of non-binary stars	*/

    int		star_circle;		/* create a circle of stars	*/
    double	cir_dist;		/* size of circle of stars	*/

    int		star_line;		/* create stars in a line	*/
    int		rnd_spacing;		/* random spacing between stars	*/

    int		do_bounce;		/* bounce the stars		*/
    int		no_speed;		/* start with no initial speed	*/
    int		few_stars;		/* start with only a few stars	*/
    int		drift;			/* let stars drift across screen*/
    int		min_angular_mom;	/* try to minimize angular mom	*/

    double	energy_fact;		/* amount of tot energy to give stars*/
    double	calc_energy_fact;	/* should this be variable	*/

    int		min_stars;		/* min # of stars before reset	*/
    int		num_add;		/* number of stars to add	*/
    int		live_stars;		/* number of live stars		*/
} sys_param_type;


/*
 * display point information
 *
 * change hash_table_bits to control the amount of memory used to record
 * displayed points.  Increasing it past 16 bits requires changing the
 * hash index table to be int's instead of shorts.
 */

#ifndef BITS
#define BITS(type)	(8 * (int)sizeof(type))
#endif


#ifndef HASH_TABLE_BITS
#define HASH_TABLE_BITS	15
#endif
#define HASH_TABLE_SIZE	(1L << (HASH_TABLE_BITS) )

#define MAX_DISP_PT	( 3L << (HASH_TABLE_BITS-2) )
#define DEFAULT_DISP_PT	( 1L << (HASH_TABLE_BITS-1) )

#define HASH_PRIME1  3238352189UL
#define HASH_PRIME2  1093769431UL
#define PT_HASH(x,y)	( ((ulong_t)(x)*HASH_PRIME1 + (ulong_t)(y)*HASH_PRIME2) >> (BITS(long) - HASH_TABLE_BITS) )


#define HASH_UNUSED	(MAX_DISP_PT+1)	/* hash entry is usable		*/
#define HASH_SEARCH	(MAX_DISP_PT+2)	/* hash entry is usable but continue */
                                        /* searching */

#define NEXT_H(a)	( ((a) >= HASH_TABLE_SIZE - 1) ? 0 : (a) + 1 )
#define PREV_H(a)	( ((a) == 0) ? HASH_TABLE_SIZE - 1 : (a) - 1 )
#define SUM_H(a,b)	( ((a) + (b) < HASH_TABLE_SIZE) ? (a) + (b) : (a) + (b) - HASH_TABLE_SIZE )
#define DIFF_H(a,b)	( ((a) >= (b)) ? (a) - (b) : (a) + HASH_TABLE_SIZE - (b) )


/*
 * display point structure.
 *
 * hstep must be large enough to hold the difference between the current
 * step number (sd.hstep) and the oldest point on the screen plus a fudge
 * factor to keep reset_hstep() from being called too often.
 *
 * the size of hstep could be made smaller by using a smaller hstep_scale
 * value.  Making it too small, however, will make the updates jerky.
 */

typedef struct
{
    XPoint	pt;
    uint_t	star: 8;		/* star that created this pt	*/
    uint_t	color: 6;		/* color in rainbow mode	*/
    uint_t	hstep: 18;		/* step number of this point	*/
#define MAX_HSTEP	((1 << 18)-1)
    int		:0;			/* no bits available		*/
} disp_point_type;
#define DISP_PT_UNUSED	(MAX_NUM_STARS)	/* unused entry flag		*/

#define NEXT(a)		( ((a) >= num_disp_pt - 1) ? 0 : (a) + 1 )
#define PREV(a)		( ((a) == 0) ? num_disp_pt - 1 : (a) - 1 )
#define SUM(a,b)	( ((a) + (b) < num_disp_pt) ? (a) + (b) : (a) + (b) - num_disp_pt )
#define DIFF(a,b)	( ((a) >= (b)) ? (a) - (b) : (a) + num_disp_pt - (b) )



/*
 * everything you need to know about the current state of the displayed
 * stars.
 */

typedef struct
{
    time_t	tstamp;			/* time since last change	*/
    int		live_erase;		/* num of stars at last erase	*/

    int		num_disp_skipped;	/* steps since last display	*/
    int		num_poll_skipped;	/* steps since last X event poll */
    int		num_visible;		/* number of visible stars	*/
    int		num_seen;		/* number of stars on screen	*/
    int		buffer_factor;		/* X buffering factor		*/

    XPoint	*last_disp;		/* last place star was displayed*/

    int		color_number;		/* what color to draw stars in	*/
    int		*star_color;		/* color of each star		*/
    int		rnd_colors[ NUM_COLORS ]; /* list of random colors	*/ 
    
    /* full screen erasing trails */
    XPoint	*points;		/* points to be displayed	*/
    XPoint	*tmp_pts;		/* points to be displayed	*/
    int		*pixels;		/* color in multi-color mode	*/
    int		max_points;		/* number of points in buffer	*/
    int		points_used;		/* points used in buffer	*/

    /* self erasing trails */
    disp_point_type	*disp_pts;	/* array of displayed points	*/
    ushort_t	*hash_index;		/* hashed table of disp_pt idx's */

    int		next_free;
    int		next_disp;
    int		next_erase;
    int		erase_hstep;
    int		updt_hstep;

    int		num_disp_pt_32;
    int		num_disp_pt_16;
    int		num_disp_pt_8;
    int		num_disp_pt_4;

    int		redraw_min_x, redraw_max_x, redraw_min_y, redraw_max_y;
    

    /* stats on what has happened */
    int		points_plotted;		/* points plotted since last erase */
    int		total_points;		/* total points plotted		*/
    double	hstep_scale;		/* scaling factor for steps	*/
    int		raw_num_steps;		/* number of steps have been done*/
    int		num_steps;		/* scaled number of steps	*/
    int		raw_hsteps;		/* wrapped number of steps	*/
    int		hsteps;			/* scaled number of steps	*/
    int		erase_disps;		/* num of disps since last erase*/

} star_disp_type;

    

#endif
