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
 * what every movement routine needs to do upon initialization.
 */

void default_init( int n, double *m, point_2d *x, point_2d *v )
{
    if( verbose > 3 )
    {
	prt_sys_const( n, x, m, v );
	prt_sys_const_delta( orig_sc, n, x, m, v );
    }
    if( verbose > 0 )
	calc_sys_const( &orig_sc, n, x, m, v );
}
