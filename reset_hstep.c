#include "xstar.h"
#include "xstar_ext.h"

/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995  Wayne Schlitt (wayne@cse.unl.edu)
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
*/


/*
 * reset the hstep values to keep from overflowing the small hstep variable
 * in star_disp_type.  
 */

void reset_hstep( star_disp_type *sd )
{
    int		i;
    int		min_hstep;
    int		hstep_used;


    if( !num_disp_pt )
	return;
    
    /*
     * find the minimum hstep
     */
    min_hstep = sd->hsteps;

    for( i = sd->next_erase; i != sd->next_free; i = NEXT(i) )
    {
	if( sd->disp_pts[i].star == DISP_PT_UNUSED )
	    continue;
	
	if( sd->disp_pts[i].hstep < min_hstep )
	    min_hstep = sd->disp_pts[i].hstep;

    }

    if( sd->erase_hstep < min_hstep )
	min_hstep = sd->erase_hstep;


    /* make sure we won't have to call reset_hstep for at least 1k steps */
    hstep_used = sd->hsteps - sd->erase_hstep;
    if( hstep_used > MAX_HSTEP - 2048 )
	sd->erase_hstep += hstep_used - (MAX_HSTEP - 2048);

    /*
     * reset the hstep variables
     */
    sd->hsteps -= min_hstep;
    sd->raw_hsteps -= min_hstep / sd->hstep_scale;
    sd->erase_hstep -= min_hstep;

    for( i = 0; i < num_disp_pt; i++ )
    {
	if( sd->disp_pts[i].star == DISP_PT_UNUSED )
	    continue;
	
	sd->disp_pts[i].hstep -= min_hstep;
    }
}


