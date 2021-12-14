/*
** XStar -- a animated n-body solver for X windows
** Copyright (C) 1995-1996  Wayne Schlitt (wayne@midwestcs.com) and others
** Version 2.2.0  12/21/96
**
** Insofar as I have any claim to this program (and I haven't much since it's
** a hack of somebody elses) the GNU General Public License would apply.
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program; if not, write to the Free Software
**  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
**
**
** XStar was based on XGrav 1.0.1 (ca 1/20/95)
** written by David Flater (dave@case50.ncsl.nist.gov)
** The README file for XGrav is in README.xgrav.
**
** XGrav was based on XSwarm 2.3 (ca 1/4/91)
** written by Jeff Butterworth (butterwo@cs.unc.edu)
**
** XSwarm was based on psychoII (date and author unknown)
** and on meltdown (date and author also unknown)
**
** psychoII was possibly based on ico (date and author unknown)
**
** ico was probably a blatant rip-off of some other program :->
*/

#include "xstar.h"
#include "xstar_ext.h"

#include <signal.h>
#include <assert.h>

/* This makes it work with virtual root window managers. */
#include "vroot.h"


#include "patchlevel.h"

/* Who needs the this header file?? */
#if defined(NEEDS_GETOPT_H)
#include <getopt.h>
#endif



int
main(argc, argv)
int argc;
char **argv;
{
    char        *geometry = NULL;


    Initialize();

    Parse_Arguments(argc, argv, &geometry);

    /* The program acts as a screen saver if timeout is non-zero. */
    if( !timeout )
    {
        Create_Window(geometry);
        Animate();
    }
    else
    {
        Change_Screen_Saver(False);
        while(TRUE)
        {
            Wait_For_Idleness();
            Create_Big_Window();
            Animate();
        }
    }

    return 0;
}


void
Quit( int signal )
{
    if( timeout )
    {
        fprintf(stderr, "Restoring screen saver state.\n");
        Change_Screen_Saver(True);
    }

    fprintf(stderr, "Terminating because of signal %d. ", signal);
    switch( signal )
    {
        case SIGHUP:    fprintf(stderr, "(Hangup)\n"); break;
        case SIGINT:    fprintf(stderr, "(Interrupt)\n"); break;
        case SIGQUIT:   fprintf(stderr, "(Quit)\n"); break;
        case SIGILL:    fprintf(stderr, "(Illegal Instruction)\n"); break;
        case SIGFPE:    fprintf(stderr, "(Floating Point Exception)\n"); break;
        case SIGKILL:   fprintf(stderr, "(Kill)\n"); break;
        case SIGBUS:    fprintf(stderr, "(Bus Error)\n"); break;
        case SIGSEGV:   fprintf(stderr, "(Segmentation Violation)\n"); break;
        case SIGTERM:   fprintf(stderr, "(Termination from Kill)\n"); break;
        case SIGSTOP:   fprintf(stderr, "(Stop not from tty)\n"); break;
        case SIGTSTP:   fprintf(stderr, "(Stop from tty)\n"); break;
#ifdef SIGXCPU
        case SIGXCPU:   fprintf(stderr, "(Exceeded CPU Time Limit)\n"); break;
#endif
        default:        fprintf(stderr, "(Unexpected Signal)\n");
    }

    exit(-signal);
}


void
Initialize()
{
    double	rand_test;
    int		i;
    int		found_large_rand;
    
    /* Get the random number generator ready. */
    srand48(time(NULL));
    sleep(1);				/* let other xstars seed w/same time*/

    for( i = 0; i < 10; i++ )
    {
	rand_test = drand48();
	if( rand_test > 1 || rand_test < 0 )
	{
	    fprintf( stderr, "Error:  drand48 returned %g, which is <0 or >1\n", rand_test );
	    fprintf( stderr, "Most likely, drand48() was not given a prototype when compiled." );
	    fprintf( stderr, "Try defining RAND48_PROTO and recompling." );
#ifdef USE_RANDOM
	    fprintf( stderr, "It is also possible that random() is not returning 31 bit integers?\n" );
#endif	
	    exit( 1 );
	}
    }

    found_large_rand = FALSE;
    for( i = 0; i < 100; i++ )
    {
	rand_test = drand48();
	if( rand_test > .5 )
	{
	    found_large_rand = TRUE;
	    break;
	}
    }

    if( !found_large_rand )
    {
	fprintf( stderr, "Error:  random number generator not working correctly.\n" );
	fprintf( stderr, "No reasonable size numbers returned from drand48().\n" );
#ifdef USE_RANDOM
	fprintf( stderr, "Possibly random() is not returning 31 bit integers?\n" );
#endif	
	exit( 1 );
    }


    signal(SIGHUP, Quit);
    signal(SIGINT, Quit);
    signal(SIGQUIT, Quit);
    signal(SIGILL, Quit);
    signal(SIGFPE, Quit);
    signal(SIGKILL, Quit);
    signal(SIGBUS, Quit);
    signal(SIGSEGV, Quit);
    signal(SIGTERM, Quit);
    signal(SIGSTOP, Quit);
    if( timeout ) signal(SIGTSTP, Quit);
#ifdef SIGXCPU
    signal(SIGXCPU, Quit);
#endif
}


/* Parse_Arguments()                                                    */
/* Analyze command line arguments and set variables accordingly.        */

void
Parse_Arguments(argc, argv, geometry)
int     argc;
char    **argv;
char    **geometry;
{
    /* These variables are used to analyze the command line parameters. */
    int         option;
    extern char *optarg;
    double	argf;
    double	argi;
    double	buffer_scale = 1;

    /* Check the command line. */
    while( (option = getopt(argc,argv,"hrt:d:g:D:b:c:C:vRMa:m:B:l:T:") )
            != EOF)
    {
        switch( option )
        {
            case 'r':
                root = TRUE;
                break;

            case 't':
                timeout = atoi(optarg);
                break;

            case 'g':
                *geometry = optarg;
                break;

            case 'd':
                display.dname = optarg;
                break;

            case 'b':
                stars = atoi(optarg);
		if( max_stars < stars + MAX_COLLAPSAR )
		    max_stars = stars + MAX_COLLAPSAR;
                if( stars < 2 )
		{
		    Usage(*argv);
		    HandleError("You need at least two stars.", FATAL);
		}
		if( max_stars > MAX_NUM_STARS )
		{
		    Usage(*argv);
		    HandleError("Too many stars specified.", FATAL );
		}
                break;

            case 'D':
                delay = atoi(optarg) * 1000; /* convert to microseconds */
                break;

            case 'c':
                star_color = optarg;
		if( rotate_colors )
		    rotate_colors = 0;
		if( multi_colors )
		    multi_colors = 0;
                break;

            case 'C':
                bg_color = optarg;
		if( rotate_colors )
		    rotate_colors = 0;
		if( multi_colors )
		    multi_colors = 0;
                break;

	    case 'v':
		verbose++;
		break;

	    case 'R':
		rotate_colors++;
		if( multi_colors )
		    multi_colors = 0;
		break;

	    case 'M':
		multi_colors++;
		if( rotate_colors )
		    rotate_colors = 0;
		break;

	    case 'a':
		argf = atof(optarg);
		if( argf <= 0 )
		{
		    Usage(*argv);
		    HandleError("-a value must be greater than zero.", FATAL );
		}
		fv *= argf;
		fv_inv = 1/fv;
		fm /= argf*argf;
		collapsar = fm*COLLAPSAR;
		break;

	    case 'm':
		if( strcmp( optarg, "euler1" ) == 0 )
		{
		    move_fna = move_stars_euler1;
		    restart_fna = restart_euler1;
		    init_fna = init_euler1;
		    move_cost = .9;
		}
		else if( strcmp( optarg, "taylor3" ) == 0 )
		{
		    move_fna = move_stars_taylor3;
		    restart_fna = restart_taylor3;
		    init_fna = init_taylor3;
		    move_cost = 1;
		}
		else if( strcmp( optarg, "ab4" ) == 0 )
		{
		    move_fna = move_stars_ab4;
		    restart_fna = restart_ab4;
		    init_fna = init_ab4;
		    move_cost = 1;
		}
		else if( strcmp( optarg, "rk4" ) == 0 )
		{
		    move_fna = move_stars_rk4;
		    restart_fna = restart_rk4;
		    init_fna = init_rk4;
		    move_cost = 4;
		}
		else if( strcmp( optarg, "gpemce8" ) == 0 )
		{
		    move_fna = move_stars_gpemce8;
		    restart_fna = restart_gpemce8;
		    init_fna = init_gpemce8;
		    move_cost = 16;
		}
		else if( strcmp( optarg, "ab7" ) == 0 )
		{
		    move_fna = move_stars_ab7;
		    restart_fna = restart_ab7;
		    init_fna = init_ab7;
		    move_cost = 1.2;
		}
		else if( strcmp( optarg, "am7" ) == 0 )
		{
		    move_fna = move_stars_am7;
		    restart_fna = restart_am7;
		    init_fna = init_am7;
		    move_cost = 2.4;
		}
		else if( strcmp( optarg, "rk4b" ) == 0 )
		{
		    move_fna = move_stars_rk4b;
		    restart_fna = restart_rk4b;
		    init_fna = init_rk4b;
		    move_cost = 4;
		}
		else
		{
		    Usage(*argv);
		    HandleError("-m value must specify a movement method.",
				FATAL );
		}
		break;

	    case 'B':
		argf = atof(optarg);
		if( argf <= 0 )
		{
		    Usage(*argv);
		    HandleError("-B value must be greater than zero.", FATAL );
		}
		buffer_scale = argf;
		break;

	    case 'l':
		argf = atof(optarg);
		if( argf <= 0 )
		{
		    Usage(*argv);
		    HandleError("-l value must be greater than zero.", FATAL );
		}
		collide_dist = argf * argf;
		break;

	    case 'T':
		argi = atoi(optarg);
		if( argi < 0 )
		{
		    Usage(*argv);
		    HandleError("-T value must be greater than or equal to zero.", FATAL );
		}
		if( argi >= MAX_DISP_PT )
		{
		    Usage(*argv);
		    HandleError("-T value must is too large.", FATAL );
		}
		num_disp_pt = argi;
		break;

		
            case 'h':
                Usage(*argv);
		exit( 0 );
		break;

	    case '?':
                Usage(*argv);
                HandleError("The command line parameters were incorrect.",
                    FATAL);
                break;
        }
    }


    /*
     * calculate the X buffering value
     *
     * this is very approximate...  It really depends on the number
     * of FLOPS that the CPU cad do and such...  The assumption is that
     * the more stars that the user selects to have, the faster the
     * computer must be, so the more buffering we can do...
     */
    if( num_disp_pt )
	max_disp_skip = (stars*stars * buffer_scale) / (2.2*move_cost*fv);
    else
	max_disp_skip = (stars*stars * buffer_scale) / (1.9*move_cost*fv);

    if( max_disp_skip > 32*32 )
	max_disp_skip = 32*32;

    if( max_disp_skip < 2*2 )
	max_disp_skip = 2*2;



    /* make sure that the # of disp pts is large enough so that
     * if they add a *bunch* of stars, we still don't have to worry
     */
    if( num_disp_pt < (MAX_NUM_STARS + 1) * 2 )
	num_disp_pt = 0;
    

    /* The screen saver is incompatible with the */
    /* root window option. */
    if( timeout > 0 )
    {
        if( root )
        {
            printf("The root window option and the screen saver option are\n");
            printf("incompatible.  The screen saver takes precedence, so\n");
            printf("that's what you're getting.\n");
        }
        root = FALSE;
    }

    /* Open the display. */
    if( !(display.dpy = XOpenDisplay(display.dname)) )
    {
        HandleError("Cannot open display.\n", FATAL);
        exit(-1);
    }

    /* Record the screen number and root window. */
    display.screen = DefaultScreen(display.dpy);
    display.root = RootWindow(display.dpy, display.screen);

    /* Set the colors. */
    display.cmap = XDefaultColormap(display.dpy, display.screen);
    if( !display.cmap ) HandleError("There was no default colormap!", FATAL);

    if( star_color == NULL )
        display.star = WhitePixel(display.dpy, display.screen);
    else
        display.star = GetColor(&display, star_color, NULL);

    if( bg_color == NULL )
        display.bg = BlackPixel(display.dpy, display.screen);
    else
        display.bg = GetColor(&display, bg_color, &(display.bg_xcolor));

    /* Set up window size. */
    if( !root && !timeout )
    {
        winW = WINWIDTH;
        winH = WINHEIGHT;
        winX = (DisplayWidth(display.dpy, display.screen)
                - winW) / 2;
        winY = (DisplayHeight(display.dpy, display.screen)
                - winH) / 2;

        if( *geometry != NULL )
	{
	    int result;

            result = XParseGeometry(*geometry, &winX, &winY, &winW, &winH);
	    if( (result & XValue) && (result & XNegative) )
		winX += DisplayWidth(display.dpy, display.screen) - winW;
	    if( (result & YValue) && (result & YNegative) )
		winY += DisplayHeight(display.dpy, display.screen) - winH;

	    if( winW < WINMIN )
		winW = WINMIN;
	    if( winH < WINMIN )
		winH = WINMIN;
	}

	root_xoff = 0;
	root_yoff = 0;
    }
    else
    {
        winX = 0;
        winY = 0;
        if( *geometry != NULL )
            XParseGeometry(*geometry, &winX, &winY, &winW, &winH);
	root_xoff = winX;
	root_yoff = winY;
	
        winW = DisplayWidth(display.dpy, display.screen);
        winH = DisplayHeight(display.dpy, display.screen);
        winX = 0;
        winY = 0;
    }
    far_dist = (winW + winH)*2.5 + 1500;
    half_winW = winW/2;
    half_winH = winH/2;
    center_x = half_winW + root_xoff;
    center_y = half_winH + root_yoff;
    min_x = 0 - center_x - 20;
    max_x = winW - center_x + 20;
    min_y = 0 - center_y - 20;
    max_y = winH - center_y + 20;

    /* If screen saving is on, then watch for events everywhere. */
    /* Traverse the window tree. */
    if( timeout ) Traverse_Tree(display.dpy, display.root);
}


/* Change_Screen_Saver()                                                */
/* Turn the server's screen saver on or off.                            */
/* This routine should be called with on=False before it is called      */
/* with on=True.                                                        */

void
Change_Screen_Saver(on)
int     on;     /* True or False */
{
    static int  timeout, interval, blanking, exposures;
    static int  set_yet = FALSE;

    if( on )
    {
        if( set_yet )
        {
            /* Restore the old settings. */
            XSetScreenSaver(display.dpy, timeout, interval, blanking,exposures);
            XActivateScreenSaver(display.dpy);
            XResetScreenSaver(display.dpy);
	    XFlush( display.dpy );
        }
    }
    else
    {
        /* Save the old settings and turn off the server's screen saver. */
        XGetScreenSaver(display.dpy, &timeout, &interval, &blanking,&exposures);
        XSetScreenSaver(display.dpy, 0, 0, DefaultBlanking, DefaultExposures);
        XResetScreenSaver(display.dpy);
        set_yet = TRUE;
    }
}


/* Traverse_Tree()                                                      */
/* Select some events from every single window that is a descendent      */
/* of "current".                                                        */

void
Traverse_Tree(display, current)
Display *display;
Window  current;
{
    Window              my_root, parent, *children;
    unsigned int        num_children;
    int         i;


    /* Watch for signs of life from the user in this window. */
    XSelectInput(display, current, ALIVE_MASK);

    /* Who are my kids? */
    XQueryTree(display, current, &my_root, &parent, &children,
        &num_children);

    /* Examine all of the children too. */
    for( i = 0; i < num_children; i++ )
        Traverse_Tree( display, children[i] );

    /* Let's not waste any memory. */
    if( num_children ) XFree( (char *) children );
}


/* Wait_For_Idleness()                                                  */
/* Wait for "timeout" seconds of user inactivity.                       */

void
Wait_For_Idleness()
{
    int         watching = TRUE;
    int         found;
    int         timer = 0;
    XEvent      event;

    while( watching )
    {
        if( timer >= timeout ) return;

        sleep(1);
        timer++;

        found = XCheckIfEvent(display.dpy, &event, Dummy_Predicate, NULL);
        if( found ) timer = 0;

        /* Flush events. */
        while( found )
        {
            switch( event.type )
            {
                /* Watch for events in new windows too. */
                case CreateNotify:
                    {
                        XCreateWindowEvent *ev = (XCreateWindowEvent *) &event;

                        /* What an insidious bug! */
                        /* I had always assumed that when I received a */
                        /* creation notification that the window in */
                        /* question would still exist.  However, the */
                        /* big screen saver window can sometimes be */
                        /* destroyed by the time we get to this point. */
                        /* We certainly don't want to do anything here */
                        /* with a window that no longer exists. */
                        if( display.win != ev->window )
                        {
                            XSelectInput(display.dpy, ev->window, ALIVE_MASK);
                            Traverse_Tree(display.dpy, ev->window);
                        }
                    }
                    break;
                default:
                    break;
            }
            /* Check for the existence of more events. */
            found = XCheckIfEvent(display.dpy, &event, Dummy_Predicate, NULL);
        }
    }
}


Bool
Dummy_Predicate(display, event, arg)
Display *display;
XEvent  *event;
char    *arg;
{
    return(True);
}


void
Create_Big_Window()
{
    /* Window Attributes */
    unsigned long               valuemask;
    XSetWindowAttributes        xswa;

    /* First time. */
    static int          first = TRUE;

    /* Cursor Stuff */
    static Cursor       cursor;

    /* Turn the cursor the color of the background. (Invisible) */
    if( first )
    {
        /* I don't care which cursor I get. */
        cursor = XCreateFontCursor(display.dpy, 0);

        XRecolorCursor(display.dpy, cursor, &(display.bg_xcolor),
                        &(display.bg_xcolor));
        first = FALSE;
    }

    /* Create a screen sized window. */
    xswa.cursor = cursor;
    xswa.background_pixel = display.bg;
    xswa.override_redirect = True;
    xswa.do_not_propagate_mask = KeyPressMask | KeyReleaseMask |
        ButtonPressMask | ButtonReleaseMask;
    valuemask = CWCursor | CWBackPixel | CWOverrideRedirect | CWDontPropagate;
    display.win = XCreateWindow(display.dpy, display.root,
        0, 0,
        winW, winH, 0,
        DefaultDepth(display.dpy, display.screen),
        InputOutput, DefaultVisual(display.dpy, display.screen),
        valuemask, &xswa);

    XMapWindow(display.dpy, display.win);

    /* Event Mask */
    XSelectInput(display.dpy, display.win, KeyPressMask);

    /* Set up the stars' graphics context. */
    display.star_gc = XCreateGC(display.dpy, display.win, 0, NULL);
    XSetForeground(display.dpy, display.star_gc, display.star);
    XSetBackground(display.dpy, display.star_gc, display.bg);

    /* Set up an erasing graphics context. */
    display.erase_gc = XCreateGC(display.dpy, display.win, 0, NULL);
    XCopyGC(display.dpy, display.star_gc, 0xffffffff, display.erase_gc);
    XSetForeground(display.dpy, display.erase_gc, display.bg);

    /* Clear the background. */
    XFillRectangle(display.dpy, display.win, display.erase_gc,
            0,0, winW, winH);
}



/* Create_Window()                                                      */
/* Create the window, making sure to set it up correctly.               */

void
Create_Window(geometry)
char    *geometry;
{
    XSetWindowAttributes xswa;
    XSizeHints          sizehint;
    XWMHints            wmhints;
    char	        wname[256];     /* Window Name */

    if( !root )
    {
        xswa.event_mask = 0;
        xswa.background_pixel = display.bg;
        display.win = XCreateWindow(display.dpy, display.root,
            winX, winY,
            winW, winH, 0,
            DefaultDepth(display.dpy, display.screen),
            InputOutput, DefaultVisual(display.dpy, display.screen),
            CWEventMask | CWBackPixel , &xswa);

        sizehint.x = winX;
        sizehint.y = winY;
        sizehint.width = winW;
        sizehint.height = winH;
        sizehint.min_width = WINMIN;
        sizehint.min_height = WINMIN;
        if( geometry != NULL )
            sizehint.flags = USPosition | USSize | PMinSize;
        else
            sizehint.flags = PPosition | PSize | PMinSize;
        XSetNormalHints(display.dpy, display.win, &sizehint);

        display.protocol_atom = XInternAtom(display.dpy, "WM_PROTOCOLS", False);
        display.kill_atom = XInternAtom(display.dpy, "WM_DELETE_WINDOW", False);
#ifdef X11R4
        XSetWMProtocols(display.dpy, display.win, &display.kill_atom, 1);
#endif

        /* Title */
        if( PATCHLEVEL )
          sprintf( wname, "XStar %s.%d", VERSION, PATCHLEVEL);
        else
          sprintf( wname, "XStar %s", VERSION);
        XChangeProperty(display.dpy, display.win,
            XA_WM_NAME, XA_STRING, 8, PropModeReplace, (uchar_t	*) wname,
            strlen(wname));

        /* Window Manager Hints (This is supposed to make input work.) */
        wmhints.flags = InputHint;
        wmhints.input = True;
        XSetWMHints(display.dpy, display.win, &wmhints);

        XMapWindow(display.dpy, display.win);
    }
    else
    {
        display.win = display.root;
	xswa.backing_store = Always;
	XChangeWindowAttributes(display.dpy, display.win,
				CWBackingStore, &xswa);
    }

    /* Event Mask */
    if( root )
        XSelectInput(display.dpy, display.win,
                        KeyPressMask | StructureNotifyMask | ExposureMask);
    else
        XSelectInput(display.dpy, display.win,
		     KeyPressMask | ButtonPressMask | StructureNotifyMask
		     | ExposureMask );

    /* Set up the stars' graphics context. */
    display.star_gc = XCreateGC(display.dpy, display.win, 0, NULL);
    XSetForeground(display.dpy, display.star_gc, display.star);
    XSetBackground(display.dpy, display.star_gc, display.bg);

    /* Set up an erasing graphics context. */
    display.erase_gc = XCreateGC(display.dpy, display.win, 0, NULL);
    XCopyGC(display.dpy, display.star_gc, 0xffffffff, display.erase_gc);
    XSetForeground(display.dpy, display.erase_gc, display.bg);

    /* Clear the background. */
    XSetWindowBackground(display.dpy, display.win, display.bg);
    XFillRectangle(display.dpy, display.win, display.erase_gc,
      0,0, winW, winH);
}



/*
** HandleEvent()
**
** process X events
*/

void
HandleEvent(event, sd)
XEvent  *event;
star_disp_type *sd;
{
    /* If the screen saver is on, then watch for signs of activity. */
    if( ((event->type == KeyPress) || (event->type == MotionNotify))
       && (timeout)
       )
	stop = TRUE;

    switch( event->type )
    {
        case ClientMessage: /* sent by f.delete from twm */
            {
                XClientMessageEvent     *ev = (XClientMessageEvent *) event;

                printf("Client message received.\n");
                if( ev->message_type == display.protocol_atom
		   && ev->data.l[0] == display.kill_atom )
                    stop = TRUE;
            }
            break;

        case ConfigureNotify:
            {
                XConfigureEvent *ev = (XConfigureEvent *) event;

                /* The screen saver should ignore resizes. */
                if( !timeout )
                {
                    winW = ev->width;
                    winH = ev->height;
		    far_dist = (winW + winH)*2.5 + 1500;
		    half_winW = winW/2;
		    half_winH = winH/2;
		    center_x = half_winW + root_xoff;
		    center_y = half_winH + root_yoff;
		    min_x = 0 - center_x - 20;
		    max_x = winW - center_x + 20;
		    min_y = 0 - center_y - 20;
		    max_y = winH - center_y + 20;
                    winX = ev->x;
                    winY = ev->y;
                }
            };
            break;

        case KeyPress:
            {
                XKeyEvent *key_event = (XKeyEvent *) event;
                char buf[128] = "";
                KeySym ks;
                XComposeStatus status;

                XLookupString(key_event,buf,128,&ks,&status);
                if( buf[0]=='q' || buf[0]=='Q' )
                    stop = TRUE;
                if( buf[0]=='d' || buf[0]=='D' )
                    add = TRUE;
                if( buf[0]=='+' )
                    add_star = TRUE;
                if( buf[0]=='-' )
                    del_star = TRUE;
                if( buf[0]=='e' || buf[0]=='E' )
                    erase = TRUE;
                if( buf[0]=='n' || buf[0]=='N' )
                    new_start = TRUE;
                if( buf[0]=='p' || buf[0]=='P' )
		    pause_updt = TRUE;
		if( buf[0]=='m' || buf[0]=='M' )
		{
		    multi_colors = !multi_colors;
		    if( rotate_colors && multi_colors )
		    {
			rotate_colors = 0;
		    }
		    if( init_colors_done == 0 )
			init_colors( colors, color_gcs, sd, NUM_COLORS );

		    do_redraw = TRUE;
		}
		if( buf[0]=='r' || buf[0]=='R' )
		{
		    rotate_colors = !rotate_colors;
		    if( rotate_colors && multi_colors )
		    {
			multi_colors = 0;
		    }
		    if( init_colors_done == 0 )
			init_colors( colors, color_gcs, sd, NUM_COLORS );

		    do_redraw = TRUE;
		}
		if( buf[0]=='\f' )
		    do_redraw = TRUE;
            }
            break;

	case Expose:
            {
		XExposeEvent	*ev = (XExposeEvent *) event;

		sd->redraw_min_x = ev->x;
		sd->redraw_max_x = ev->x + ev->width;
		sd->redraw_min_y = ev->y;
		sd->redraw_max_y = ev->y + ev->height;

		screen_exposed = TRUE;
	    }
	    break;
	    
        default:
            break;
    }
}


#ifdef USE_USLEEP

static void alarmhandler()
{
}

void sleepms(msec)
int msec;
{
    struct itimerval value,ovalue;
    struct sigvec vec;
    long savemask, sigblock(), sigpause();

    vec.sv_handler = alarmhandler;
    vec.sv_mask = 0x0;
    vec.sv_flags = 0;
    sigvector(SIGALRM, &vec, &vec); /* Set up alarmhandler for SIGALRM */
    savemask = sigblock((long)(1L << (SIGALRM - 1)));

    value.it_interval.tv_sec = 0;
    value.it_interval.tv_usec = 0;
    value.it_value.tv_sec = msec/1000;
    value.it_value.tv_usec = (msec%1000)*1000;
    setitimer(ITIMER_REAL,&value,&ovalue);

    (void)sigpause(0L);
    (void)sigsetmask(savemask);

    sigvector(SIGALRM, &vec, NULL); /* Restore previous signal handler */
}

void usleep(us)
unsigned us;
{
    sleepms(us / 1000);
}

#else

/*  Put the process to sleep for a while.  */
void nap(sec,usec)
long sec, usec;
{
    static struct timeval tv;

    tv.tv_sec = sec;
    tv.tv_usec = usec;
    select(0, 0, 0, 0, &tv);
}

#endif


void
Usage(program)
char *program;
{
    printf("%s [options]  where options are listed below\n", program);
    printf("-h                display this message\n");
    printf("-b stars          number of stars  (2 to %d)\n",
	   MAX_NUM_STARS - MAX_COLLAPSAR );
    printf("-r                use root window\n");
    printf("-g geometry       window geometry\n");
    printf("-d host:display   X server to connect to\n");
    printf("-t timeout        screen saved after 'timeout' seconds\n");
    printf("-D delay          delay between updates (milliseconds)\n");
    printf("-c star_clr       star color\n");
    printf("-C bg_color       background color\n");
    printf("-R                Rotate the star colors  (Default unless -c or -C is used)\n");
    printf("-M                multiple colors, one per star\n");
    printf("-a float          accuracy of position calculations.   Larger values\n" );
    printf("                  increase accuracy.  Must be greater than zero.\n" );
    printf("-m [euler1|taylor3|rk4|rk4b|gpemce8|ab7|am7]     method used to calc points\n" );
    printf("-l float          distance between stars before they collide\n");
    printf("-T num_pts        number of points to display as star trails\n");
    printf("-B buf_factor     amount of buffering to be done with the X server.\n" );
    printf("-v                enable display of debug output\n");
    printf("\nDefaults are equivalent to the following command:\n");
    printf("%s -b %d -g =%dx%d -D %d -c White -C Black -R -a 1.0 -m ab7 -l %.1f -T %ld -B 1.0\n",
	   program, STARS, WINWIDTH, WINHEIGHT, DELAY,
	   DEFAULT_COLLIDE, DEFAULT_DISP_PT );
    printf("\nCommands available while running are:\n" );
    printf("  d  create a collapsar (gravity well).\n");
    printf("  e  erase trails.\n");
    printf("  +  add a star to the system.\n");
    printf("  -  delete a star from the system.\n");
    printf("  n  get a new set of stars.\n");
    printf("  m  toggle multi-color mode.\n");
    printf("  r  toggle rainbow mode.\n");
    printf("  p  pause the updating.  Press p again to start.\n" );
    printf("  ^L redraw the screen.\n" );
    printf("  q  quit.\n\n");
}


void
HandleError(description, degree)
char    *description;
int     degree;
{
    fprintf(stderr, "An error has occurred.  The description is below...\n");
    fprintf(stderr, "%s\n", description);

    if( degree == FATAL )
    {
        fprintf(stderr, "Program aborting...\n");
        exit(-1);
    }
}

long
GetColor(display, color, final_color)
disp            *display;
char            *color;
XColor          *final_color;
{
    XColor      cdef;
    char        error_str[STD_STR];

    if( !XParseColor(display->dpy, display->cmap, color, &cdef) ||
        !XAllocColor(display->dpy, display->cmap, &cdef) )
    {
        sprintf(error_str, "Color \"%s\" wasn't found.", color);
        HandleError(error_str, FATAL);
    }

    /* Copy the final color. */
    if( final_color != NULL ) *final_color = cdef;

    return(cdef.pixel);
}


/* Check the return code of malloc(). */
void *
Safe_Malloc(bytes)
int bytes;
{
    void *pointer;

    pointer = (void *) malloc(bytes);
    if( NULL == pointer )
    {
        fprintf(stderr, "Error allocating %d bytes.\n", bytes);
        exit(-1);
    }

    return(pointer);
}

void *
Safe_Realloc(ptr, bytes)
void *ptr;
int bytes;
{
    void *pointer;

    pointer = (void *) realloc(ptr, bytes);
    if( NULL == pointer )
    {
        fprintf(stderr, "Error allocating %d bytes.\n", bytes);
        exit(-1);
    }

    return(pointer);
}
