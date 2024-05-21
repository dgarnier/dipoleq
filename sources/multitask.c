#include "multitask.h"

#if MULTITASK
#ifdef __MWERKS__
#include <MacTypes.h>
//#include <Quickdraw.h>
#include <Events.h>
//#include <Resources.h>
#include <SIOUX.h>
#include <Appearance.h>

/* Cursor related functions */
//Boolean InitAnimatedCursors( short acurID );
//void ReleaseAnimatedCursors( void );
//void SpinMyCursor( void );

/* Private globals */
static long			gNextCheck;			/* Used to hold tick counts */
static long			gTickSlice = 15;	/* 1/30th second in ticks (each tick = 1/60th sec) */


void doMultiTask( long sleepTime )
{
	static UInt32 spin = 0;
	extern Boolean SIOUXQuitting;
	EventRecord myEvent;
	OSStatus err;

 	if( TickCount() > gNextCheck )   					/* Time to check for events again? */
    {
     	if( WaitNextEvent( everyEvent, &myEvent, sleepTime, NULL ) )
      	{
      		SetThemeCursor( kThemeArrowCursor );
      		/* Restore arrow while we're handling a real event */
      //		SetCursor(NULL);

      		/* Need to do something with the event if we got one */
			/* Add additional event handling code here if you need it */
			SIOUXHandleOneEvent( &myEvent );
		//	if( SIOUXQuitting )
    	//		sig_die("User interrupt; execution stopped", 1);		/* Graceful quit */
 		}

    	err = SetAnimatedThemeCursor ( kThemeSpinningCursor, spin++);

    	gNextCheck = TickCount() + gTickSlice;			/* Reset the tick count */
    }
}





#endif /* __MWERKS__ */
#endif /* MULTITASK */
