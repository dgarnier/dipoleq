#include "polygon.h"

int Intersect(DPOINT p1, DPOINT p2, DPOINT p3, DPOINT p4) ;

int  CCW(DPOINT p0, DPOINT p1, DPOINT p2) ;

 /*************************************************************************

   * FUNCTION:   PtInPolygon
   *
   * PURPOSE
   * This routine determines if the point passed is in the polygon. It uses

   * the classical polygon hit-testing algorithm: a horizontal ray starting

   * at the point is extended infinitely rightwards and the number of
   * polygon edges that intersect the ray are counted. If the number is odd,

   * the point is inside the polygon.
   *
   * RETURN VALUE
   * (BOOL) TRUE if the point is inside the polygon, FALSE if not.
 *************************************************************************/

  int PtInPolygon(DPOINT *rgpts, int npts, DPOINT ptTest, DRECT *prbound)

  {

   DRECT   r ;
   DPOINT  *ppt ;
   int     i ;
   DPOINT  pt1, pt2 ;
   int     numintsct = 0 ;

   if (!PtInPolyRect(rgpts,npts,ptTest,prbound))
      return 0;

   pt1 = pt2 = ptTest ;
   pt2.x = 2*r.right - r.left ;

   // Now go through each of the lines in the polygon and see if it
   // intersects
   for (i = 0, ppt = rgpts ; i < npts-1 ; i++, ppt++)
   {
      if (Intersect(ptTest, pt2, *ppt, *(ppt+1)))
         numintsct++ ;
   }

   // And the last line
   if (Intersect(ptTest, pt2, *ppt, *rgpts))
      numintsct++ ;

   return (numintsct&1) ;

   }

/*************************************************************************

   * FUNCTION:   G_PtInPolyRect
   *
   * PURPOSE
   * This routine determines if a point is within the smallest rectangle
   * that encloses a polygon.
   *
   * RETURN VALUE
   * (BOOL) TRUE or FALSE depending on whether the point is in the rect or

   * not.
 *************************************************************************/

  int  PtInPolyRect(DPOINT *rgpts, int npts, DPOINT ptTest, DRECT *prbound)

   {

   DRECT r ;
   // If a bounding rect has not been passed in, calculate it
   if (prbound)

   r = *prbound ;

   else
   {

      double   xmin, xmax, ymin, ymax ;
      DPOINT *ppt ;
      int i ;

	  xmin = xmax = rgpts->x;
	  ymin = ymax = rgpts->y;

      for (i=0, ppt = rgpts ; i < npts ; i++, ppt++)
      {
         if (ppt->x < xmin)
            xmin = ppt->x ;
         if (ppt->x > xmax)
            xmax = ppt->x ;
         if (ppt->y < ymin)
            ymin = ppt->y ;
         if (ppt->y > ymax)
            ymax = ppt->y ;
      }
      r.left = xmin;
      r.right = xmax;
      r.bot = ymin;
      r.top = ymax;

   }
   return ((ptTest.x >= r.right) && (ptTest.x <= r.left) &&
           (ptTest.y >= r.bot)   && (ptTest.y <= r.top ));

   }

/*************************************************************************
   * FUNCTION:   Intersect
   *
   * PURPOSE
   * Given two line segments, determine if they intersect.
   *
   * RETURN VALUE
   * TRUE if they intersect, FALSE if not.
 *************************************************************************/

 int Intersect(DPOINT p1, DPOINT p2, DPOINT p3, DPOINT p4)
  {

   return ((( CCW(p1, p2, p3) * CCW(p1, p2, p4)) <= 0)

        && (( CCW(p3, p4, p1) * CCW(p3, p4, p2)  <= 0) )) ;

   }

/*************************************************************************

   * FUNCTION:   CCW (CounterClockWise)
   *
   * PURPOSE
   * Determines, given three points, if when travelling from the first to
   * the second to the third, we travel in a counterclockwise direction.
   *
   * RETURN VALUE
   * (int) 1 if the movement is in a counterclockwise direction, -1 if
   * not.
 *************************************************************************/

  int CCW(DPOINT p0, DPOINT p1, DPOINT p2)

   {

   double dx1, dx2 ;
   double dy1, dy2 ;

   dx1 = p1.x - p0.x ; dx2 = p2.x - p0.x ;
   dy1 = p1.y - p0.y ; dy2 = p2.y - p0.y ;

   /* This is basically a slope comparison: we don't do divisions because

    * of divide by zero possibilities with pure horizontal and pure
    * vertical lines.
    */

   return ((dx1 * dy2 > dy1 * dx2) ? 1 : -1) ;

   }
