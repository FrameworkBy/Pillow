#include "Python.h"
#include "Imaging.h"


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

Imaging ImagingNlmFilter(Imaging imOut, Imaging imIn, float sigma)
{
    for (y = 0; y < im->ysize; y++)
		for (x = 0; x < im->xsize; x++)
			imOut[x][y] = imIn[x][y]+sigma;

    return imOut;
}