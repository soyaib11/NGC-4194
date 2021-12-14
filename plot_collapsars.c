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



/* plot the collapsars on the screen */
void plot_collapsars( point_2d *cur, double *m, sys_param_type *sp, star_disp_type *sd, disp *display )
{

    int		i;

    int		x, y;

    int		num_pts;

    XPoint	*disp_buf;


    /* find a good scratch array */
    if( num_disp_pt )
	disp_buf = sd->points;
    else
	disp_buf = sd->tmp_pts;
    
    

    /* draw the collapsars */
    
    for( i = 0; i < max_stars; i++ )
    {
	if( m[i] < collapsar )
	    continue;

	
	num_pts = 0;
	for( x = -1; x < 2; x++ )
	    for( y = -1; y < 2; y++ )
	    {
		disp_buf[num_pts].x = cur[i].x + center_x + x;
		disp_buf[num_pts].y = cur[i].y + center_y + y;
		num_pts++;
		
	    }
	
	if( rotate_colors )
	    XDrawPoints( display->dpy, display->win,
			color_gcs[ sd->color_number ],
			disp_buf, num_pts,
			CoordModeOrigin );
	else if( multi_colors )
	    XDrawPoints( display->dpy, display->win,
			color_gcs[ sd->star_color[ i ]],
			disp_buf, num_pts,
			CoordModeOrigin );
	else 
	    XDrawPoints( display->dpy, display->win,
			display->star_gc, disp_buf, num_pts,
			CoordModeOrigin );
	
    }
}

