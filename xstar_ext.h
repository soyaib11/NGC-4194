#ifndef include_xstar_ext
#define include_xstar_ext


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com) and others
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/




#ifdef RAND48_PROTO
double drand48(void);
long lrand48(void);
void srand48( long seedval );
#endif

#ifdef USE_RANDOM
long  random(void);
void  srandom(int seed);
#endif


/*--------------------------- Function Prototypes ---------------------------*/

extern void
Quit(
#ifdef __STDC__
int signal
#endif
);

extern void    Initialize();

extern void
Parse_Arguments(
#ifdef __STDC__
int  argc,
char **argv,
char **geometry
#endif
);

extern void
Change_Screen_Saver(
#ifdef __STDC__
int
#endif
);

extern void
Traverse_Tree(
#ifdef __STDC__
Display *,
Window
#endif
);

extern void Wait_For_Idleness();

extern Bool
Dummy_Predicate(
#ifdef __STDC__
Display *,
XEvent *,
char *
#endif
);

extern void Create_Big_Window();

extern void
Create_Window(
#ifdef __STDC__
char *
#endif
);

extern void init_stars( point_2d *cur, point_2d *prev, double *m,
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd );

extern void new_star( point_2d *cur, point_2d *prev, double *m, 
		point_2d *vk[], point_2d *ak[],
		int max_stars, sys_param_type *sp, star_disp_type *sd,
		int add_collapsar );


extern void set_sys_param( sys_param_type *sp, int max_stars );
    
extern void set_xmva( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void set_va( point_2d *cur, point_2d *prev, double *m,
	    point_2d *vk[], point_2d *ak[],
	    sys_param_type *sp,
	    int i );

extern void set_star_disp( star_disp_type *sd, sys_param_type *sp );

extern void set_buffer_factor( star_disp_type *sd, sys_param_type *sp );


extern void init_sys_1a( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void init_sys_1b( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void init_sys_1c( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void init_sys_4( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void init_sys_8( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

extern void init_sys_8b( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );

void dump_sys( point_2d *cur, point_2d *prev, double *m,
	      point_2d *vk[], point_2d *ak[],
	      int max_stars, sys_param_type *sp );


extern void init_colors( XColor *colors, GC *color_gcs, star_disp_type *sd, int num_elem );

extern void restart_euler1( int n, double *m, point_2d *x, point_2d *v );

extern void init_euler1( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_euler1( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_taylor3( int n, double *m, point_2d *x, point_2d *v );

extern void init_taylor3( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_taylor3( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_ab4( int n, double *m, point_2d *x, point_2d *v );

extern void init_ab4( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_ab4( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_rk4( int n, double *m, point_2d *x, point_2d *v );

extern void init_rk4( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_rk4( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_gpemce8( int n, double *m, point_2d *x, point_2d *v );

extern void init_gpemce8( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_gpemce8( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_ab7( int n, double *m, point_2d *x, point_2d *v );

extern void init_ab7( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_ab7( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_tam5( int n, double *m, point_2d *x, point_2d *v );

extern void init_tam5( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_tam5( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_am7( int n, double *m, point_2d *x, point_2d *v );

extern void init_am7( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_am7( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_tab4( int n, double *m, point_2d *x, point_2d *v );

extern void init_tab4( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_tab4( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void restart_rk4b( int n, double *m, point_2d *x, point_2d *v );

extern void init_rk4b( int n, double *m, point_2d *x, point_2d *v );

extern void move_stars_rk4b( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void default_init( int n, double *m, point_2d *x, point_2d *v );


extern void check_bounce( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );


extern void do_collide( int i, int j, point_2d *cur, point_2d *prev, double *m,
		point_2d *vk[], point_2d *ak[],
		sys_param_type *sp, star_disp_type *sd,
		char *file, int line 
		);

extern double kinetic_energy( int i, int n, point_2d *cur, double *m, point_2d *v );

extern double binding_energy( int i, int n, point_2d *cur, double *m, point_2d *v );


extern void calc_sys_const( sys_const *sc,
			   int n, point_2d *cur, double *m,  point_2d *v );
extern void do_prt_sys_const( int n, point_2d *cur, double *m,  point_2d *v, char *file, int line );
extern void do_prt_sys_const_delta( sys_const sc0, int n, point_2d *cur, double *m,  point_2d *v, char *file, int line );

extern void update_screen( sys_param_type *sp, star_disp_type *sd, disp *display );

extern void clear_disp_pt( star_disp_type *sd );

extern void redraw_screen(  point_2d *cur, double *m, sys_param_type *sp, star_disp_type *sd, disp *display );

extern void set_redraw_full( star_disp_type *sd );

extern void plot_collapsars( point_2d *cur, double *m, sys_param_type *sp, star_disp_type *sd, disp *display );


extern void reset_hstep( star_disp_type *sd );


extern void Animate();

extern void
HandleEvent(
#ifdef __STDC__
XEvent *,
star_disp_type *sd
#endif
);

#ifdef USE_USLEEP
extern void alarmhandler();

extern void sleepms(
#ifdef __STDC__
int msec
#endif
);

extern void usleep(
#ifdef __STDC__
unsigned us
#endif
);

#else
extern void    nap();
#endif

extern void    Usage(
#ifdef __STDC__
char *program
#endif
);

extern void    HandleError(
#ifdef __STDC__
char    *description,
int     degree
#endif
);

extern long
GetColor(
#ifdef __STDC__
disp            *display,
char            *color,
XColor          *final_color
#endif
);

extern void *
Safe_Malloc(
#ifdef __STDC__
int bytes
#endif
);

extern void *
Safe_Realloc(
#ifdef __STDC__
void *ptr,
int bytes
#endif
);


/*----------------------------- Global Variables ----------------------------*/

/* X related */
extern int		winX, winY;
extern int		root_xoff, root_yoff;
extern double		center_x, center_y;
extern uint_t		winW, winH;
extern int		half_winW, half_winH;
extern int		max_x, min_x, max_y, min_y;
extern uint_t		far_dist;
extern disp		display;

/* animation related */
extern int     stars;          /* number of stars */
extern int     max_stars;    /* maximum number of stars */
extern int     max_disp_skip;    /* max num of displays to skip */
extern int     delay;   /* delay between updates, in microseconds */
extern int     timeout;            /* time in seconds before screen saving */
extern char    *star_color;
extern char    *bg_color;
extern XColor	colors[NUM_COLORS];
extern GC	color_gcs[NUM_COLORS];
extern int	verbose;
extern int	rotate_colors;
extern int	init_colors_done;
extern int	multi_colors;
extern int	num_disp_pt;

extern int     stop;
extern int     add;           /* flag to add a collapsar */
extern int     add_star;      /* flag to add a star */
extern int     del_star;      /* flag to delete a star */
extern int     erase;          /* flag to erase trails */
extern int     new_start;      /* flag to start a new system trails */
extern int	pause_updt;	/* flag to pause the updates	*/
extern int	screen_exposed;		/* flag to update the screen	*/
extern int	do_redraw;		/* flag to redraw the screen	*/


extern int     root;           /* display in root window */
extern time_t	last_init;

extern double	collapsar;
extern double	fv;
extern double	fv_inv;
extern double	fm;

extern sys_param_type	orig_sp;
extern sys_const orig_sc;
extern star_disp_type	orig_sd;

extern int num_collide;
extern int num_edge;

extern double	collide_dist;

extern void (*move_fna)( point_2d *prev, point_2d *cur, double *m, 
		  point_2d *vk[], point_2d *ak[],
		  int max_stars, sys_param_type *sp, star_disp_type *sd );

extern void (*init_fna)( int n, double *m, point_2d *x, point_2d *v );

extern void (*restart_fna)( int n, double *m, point_2d *x, point_2d *v );

extern double move_cost;

#endif
