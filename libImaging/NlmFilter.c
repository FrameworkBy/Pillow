#include "Python.h"
#include "Imaging.h"

Imaging ImagingNlmFilter(Imaging imOut, Imaging imIn, float sigma)
{
	int x, y;
	ImagingCopyInfo(imOut, imIn);
	if (imIn->image8)
		for (y = 0; y < imIn->ysize; y++)
			for (x = 0; x < imIn->xsize; x++)
				imOut->image8[x][y] = imIn->image8[x][y];

    return imOut;
}