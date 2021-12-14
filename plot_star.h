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
 * this routine is common between all the movement routines, but it's
 * too expensive to make it a function call, so it is inlined via the
 * #include command.
 */

	/* if the star is on the screen, then plot it */
	if( xpt > min_x && ypt > min_y && xpt < max_x && ypt < max_y )
	{
	    if( m_i >= collapsar )
		continue;

	    ++sd->num_visible;
	    ++sd->num_seen;
	    
	    if( xpt != sd->last_disp[i].x || ypt != sd->last_disp[i].y )
	    {
		sd->last_disp[i].x = xpt;
		sd->last_disp[i].y = ypt;
		
		xpt += center_x;
		ypt += center_y;
		    
		if( num_disp_pt )
		{
		    disp_point_type	*new_pt = &sd->disp_pts[sd->next_free];

		    int		hash_loc;
		    int		hashed_idx;
		    int		first_search = -1;


#ifdef DUMP_STAR_DATA
		    printf( "%d  %d  %.16e %.16e  %.16e %.16e  %.16e %.16e\n", i, sd->raw_num_steps, cur[i].x, cur[i].y, v_0[i].x, v_0[i].y, ak[0][i].x, ak[0][i].y );
#endif


		    new_pt->pt.x  = xpt;
		    new_pt->pt.y  = ypt;
		    new_pt->star  = i;
		    new_pt->color = sd->color_number;
		    new_pt->hstep = sd->hsteps;


		    /*
		     * check to see if there was already a star plotted here
		     */
		    for( hash_loc = PT_HASH( xpt, ypt ),
			hashed_idx = sd->hash_index[ hash_loc ];
			
			hashed_idx != HASH_UNUSED;
			
			hash_loc = NEXT_H( hash_loc ),
			hashed_idx = sd->hash_index[ hash_loc ] )
		    {
			if( hashed_idx == HASH_SEARCH )
			{
			    if( first_search == -1 )
				first_search = hash_loc;

			    continue;
			}
			
			if( sd->disp_pts[ hashed_idx ].pt.x != xpt
			   || sd->disp_pts[ hashed_idx ].pt.y != ypt
			   )
			    continue;


			sd->disp_pts[ hashed_idx ].star = DISP_PT_UNUSED;

			/*
			 * If we have found a search entry, then we should
			 * really use it instead of this entry and we should
			 * delete this entry.  But this case happens so 
			 * rarely (literally one in a million times) that
			 * it isn't worth it...  
			 */
			first_search = -1;

			break;
		    }

		    if( first_search == -1 )
			sd->hash_index[ hash_loc ] = sd->next_free;
		    else
			sd->hash_index[ first_search ] = sd->next_free;
		    
		    
		    sd->next_free = NEXT( sd->next_free );

		}
		else
		{
		    if( multi_colors )
			sd->pixels[sd->points_used] = sd->star_color[ i ];
		    
		    sd->points[sd->points_used].x = xpt;
		    sd->points[(sd->points_used)++].y = ypt;
		}
	    }
	}
	else
	{
	    /* see if the star has left the system */
	    
	    if( abs(xpt) + abs(ypt) > far_dist )
	    {
		num_edge++;
		--sp->live_stars;
		set_buffer_factor( sd, sp );
		if( verbose > 1 )
		{
		    double	k, u;
		    double	dx, dy;
		    double	tdist;
		
		    dx = cur[i].x;
		    dy = cur[i].y;
		    tdist = sqrt( dx*dx + dy*dy );
		    
		    k = kinetic_energy( i, max_stars, cur, m, v_0 );
		    u = binding_energy( i, max_stars, cur, m, v_0 );
		    
		    if( k + u < 0 )
			fprintf( stderr,
				"%2d  %6.1f,%6.1f  %7.4f,%7.4f %6.3f  %6.2f %9.6f %9.6f %9.6f\n",
				i, cur[i].x, cur[i].y, fv*v_0[i].x, fv*v_0[i].y, m[i]/fm,
				tdist, k*fv*fv/fm, u*fv*fv/fm, (k + u)*fv*fv/fm );
	    
		    if( verbose > 1 )
			fprintf( stderr, "%s:%d  live_stars=%d  star %d  |(%g,%g)|=%d\n", __FILE__, __LINE__, sp->live_stars, i, cur[i].x, cur[i].y, xpt*xpt+ypt*ypt );
		}
		
		m[i] = m_i = 0.0;
		cur[i].x = cur[i].y = prev[i].x = prev[i].y = -1.0;
		for( k = 0; k < K; k++ )
		{
		    vk[k][i].x = vk[k][i].y = 0;
		    ak[k][i].x = ak[k][i].y = 0;
		}	    

		init_fna( max_stars, m, cur, v_0 );
	    }
	}
