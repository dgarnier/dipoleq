/************************************************************************************
**
**  polygon.h
**
**    
**
************************************************************************************/

typedef struct {
	double x, y;
} DPOINT;

typedef struct {
	double left, top, right, bot;
} DRECT;

int PtInPolygon(DPOINT *rgpts, int npts, DPOINT ptTest, DRECT *prbound) ;

int PtInPolyRect(DPOINT *rgpts, int npts, DPOINT ptTest, DRECT *prbound) ;