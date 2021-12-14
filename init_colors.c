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



/* straight out of Foley and Van Dam... */
static void hsv_to_rgb( double h, double s, double v,
		double *r, double *g, double *b )
{
    double f;
    double p;
    double q;
    double t;
    int sextant;
    
    if( h == 360.0 )
	h = 0.0;

    h /= 60.0;
    sextant = h;
    f = h - sextant;
    p = v * (1.0 - s);
    q = v * (1.0 - s * f);
    t = v * (1.0 - s * (1.0 - f));
    
    switch( sextant )
    {
    case 0 :
	*r = v;
	*g = t;
	*b = p;
	break;
	
    case 1 :
	*r = q;
	*g = v;
	*b = p;
	break;
	
    case 2 :
	*r = p;
	*g = v;
	*b = t;
	break;
	
    case 3 :
	*r = p;
	*g = q;
	*b = v;
	break;
	
    case 4 :
	*r = t;
	*g = p;
	*b = v;
	break;
	
    case 5 :
	*r = v;
	*g = p;
	*b = q;
	break;
    }
}


/*
 * initialize the color array to form a rainbow
 */

void init_colors( XColor *colors, GC *color_gcs, star_disp_type *sd, int num_elem )
{
    int			c;

    double		h, s = 1, v = 1;
    double		r, g, b;

    double		h_inc = 360. / NUM_COLORS;
    
    char        error_str[STD_STR];


    for( c = 0; c < NUM_COLORS; c++ )
    {
	h = c * h_inc;
	hsv_to_rgb( h, s, v, &r, &g, &b );


	/* allocate another color */
	colors[c].red = r * MAX_CVAL;
	colors[c].green = g * MAX_CVAL;
	colors[c].blue = b * MAX_CVAL;
	if( !XAllocColor(display.dpy, display.cmap, &colors[ c ]) )
	{
	    sprintf(error_str, "Color %d=(%d,%d,%d) couldn't be allocated.",
		    c, colors[c].red, colors[c].green,
		    colors[c].blue );
	    HandleError(error_str, FATAL);
	}
#if 0
	fprintf( stderr, "%s:%d  color[%d]=(%d,%d,%d)=%08X\n", __FILE__, __LINE__, c, colors[c].red, colors[c].green, colors[c].blue, colors[c].pixel );
#endif
	
	color_gcs[c] = XCreateGC(display.dpy, display.win, 0, NULL);
	XSetForeground( display.dpy, color_gcs[c], colors[c].pixel);
    }

    init_colors_done = 1;
}

