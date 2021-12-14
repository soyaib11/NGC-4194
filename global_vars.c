#include "xstar.h"
#include "xstar_ext.h"


/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com) and others
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/



/* X related */
int		winX, winY;
int		root_xoff, root_yoff;
double		center_x, center_y;
uint_t		winW, winH;
int		half_winW, half_winH;
int		max_x, min_x, max_y, min_y;
uint_t		far_dist;
disp		display;

/* animation related */
int     stars = STARS;          /* number of stars */
int     max_stars = STARS+MAX_COLLAPSAR;    /* maximum number of stars */
int     max_disp_skip;			/* max num of displays to skip */
int     delay = DELAY * 1000;   /* delay between updates, in microseconds */
int     timeout = 0;            /* time in seconds before screen saving */
char    *star_color = NULL;
char    *bg_color = NULL;
XColor	colors[NUM_COLORS];
GC	color_gcs[NUM_COLORS];
int	verbose = 0;
int	rotate_colors = 1;
int	init_colors_done = 0;
int	multi_colors = 0;
int	num_disp_pt = DEFAULT_DISP_PT;

int     stop = FALSE;
int     add  = FALSE;           /* flag to add a collapsar */
int     add_star  = FALSE;      /* flag to add a star */
int     del_star  = FALSE;      /* flag to delete a star */
int     erase = FALSE;          /* flag to erase trails */
int     new_start = FALSE;      /* flag to start a new system trails */
int	pause_updt = FALSE;	/* flag to pause the updates	*/
int	screen_exposed;		/* flag to update the screen	*/
int	do_redraw;		/* flag to redraw the screen	*/

int     root = FALSE;           /* display in root window */
time_t	last_init = 0;

double	collapsar = (FM*COLLAPSAR);
double	fv = FV;
double	fv_inv = 1./FV;
double	fm = FM;

sys_param_type	orig_sp;
sys_const orig_sc;
star_disp_type  orig_sd;

int num_collide;
int num_edge;

double	collide_dist = DEFAULT_COLLIDE*DEFAULT_COLLIDE;

void (*move_fna)( point_2d *prev, point_2d *cur, double *m, 
		 point_2d *vk[], point_2d *ak[],
		 int max_stars, sys_param_type *sp, star_disp_type *sd )
    = move_stars_ab7;
    
void (*init_fna)( int n, double *m, point_2d *x, point_2d *v )
    = init_ab7;

void (*restart_fna)( int n, double *m, point_2d *x, point_2d *v )
    = restart_ab7;

double move_cost = 1.2;
