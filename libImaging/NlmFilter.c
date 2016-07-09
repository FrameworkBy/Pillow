#include "Python.h"
#include "Imaging.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>




#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )


#define dTiny 1e-10
#define fTiny 0.00000001f
#define fLarge 100000000.0f

///// LUT tables
#define LUTMAX 30.0
#define LUTMAXM1 29.0
#define LUTPRECISION 1000.0

void  wxFillExpLut(float *lut,int size);        // Fill exp(-x) lut

float wxSLUT(float dif,float *lut);                     // look at LUT

void fpClear(float *fpI,float fValue, int iLength);

float fiL2FloatDist1 ( float * u0, float *u1, int i0, int j0,int i1,int j1,int radius, int width0, int width1);

float fiL2FloatDist2 (float **u0, float **u1,int i0,int j0,int i1,int j1,int radius,int channels, int width0, int width1);

void nlmeans_ipol(int iDWin,                    // Half size of comparison window
                  int iDBloc,           // Half size of research window
                  float fSigma,         // Noise parameter
                  float fFiltPar,       // Filtering parameter
                  float **fpI,          // Input
                  float **fpO,          // Output
                  int iChannels, int iWidth,int iHeight);


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

void fpClear(float *fpI,float fValue, int iLength) {
  int ii;
    for (ii=0; ii < iLength; ii++) fpI[ii] = fValue;
}






// LUT tables
void  wxFillExpLut(float *lut, int size) {
  int i;
    for (i=0; i< size; i++)   lut[i]=   expf( - (float) i / LUTPRECISION);
}




float wxSLUT(float dif, float *lut) {

    if (dif >= (float) LUTMAXM1) return 0.0;

    int  x= (int) floor( (double) dif * (float) LUTPRECISION);

    float y1=lut[x];
    float y2=lut[x+1];

    return y1 + (y2-y1)*(dif*LUTPRECISION -  x);
}






float fiL2FloatDist1(float *u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width0, int width1) {

    float dist=0.0;
    int s, r;
    for (s=-radius; s<= radius; s++) {

        int l = (j0+s)*width0 + (i0-radius);
        float *ptr0 = &u0[l];

        l = (j1+s)*width1 + (i1-radius);
        float *ptr1 = &u1[l];

        for (r=-radius; r<=radius; r++,ptr0++,ptr1++) {
            float dif = (*ptr0 - *ptr1);
            dist += (dif*dif);
        }

    }

    return dist;
}




float fiL2FloatDist2(float **u0,float **u1,int i0,int j0,int i1,int j1,int radius,int channels, int width0, int width1) {

    float dif = 0.0f;


    for (int ii=0; ii < channels; ii++) {

        dif += fiL2FloatDist1(u0[ii],u1[ii],i0,j0,i1,j1,radius,width0, width1);

    }


    return dif;
}

void nlmeans_ipol(int iDWin,            // Half size of patch
                  int iDBloc,           // Half size of research window
                  float fSigma,         // Noise parameter
                  float fFiltPar,       // Filtering parameter
                  float **fpI,          // Input
                  float **fpO,          // Output
                  int iChannels, int iWidth,int iHeight) {




    // length of each channel
    int iwxh = iWidth * iHeight;


    //  length of comparison window
    int ihwl = (2*iDWin+1);
    int iwl = (2*iDWin+1) * (2*iDWin+1);
    int icwl = iChannels * iwl;



    // filtering parameter
    float fSigma2 = fSigma * fSigma;
    float fH = fFiltPar * fSigma;
    float fH2 = fH * fH;

    // multiply by size of patch, since distances are not normalized
    fH2 *= (float) icwl;



    // tabulate exp(-x), faster than using directly function expf
    int iLutLength = (int) rintf((float) LUTMAX * (float) LUTPRECISION);
    float *fpLut = (float*)malloc(sizeof(float)*iLutLength);
    wxFillExpLut(fpLut,iLutLength);




    // auxiliary variable
    // number of denoised values per pixel
    float *fpCount = (float*)malloc( sizeof(float)*iwxh);
    fpClear(fpCount, 0.0f,iwxh);




    // clear output
    for (int ii=0; ii < iChannels; ii++) fpClear(fpO[ii], 0.0f, iwxh);




    // PROCESS STARTS
    // for each pixel (x,y)
//#pragma omp parallel shared(fpI, fpO)
    {


//#pragma omp for schedule(dynamic) nowait

        for (int y=0; y < iHeight ; y++) {


            // auxiliary variable
            // denoised patch centered at a certain pixel
            float **fpODenoised = (float**)malloc( sizeof(float*)*iChannels);
            for (int ii=0; ii < iChannels; ii++) fpODenoised[ii] = (float*)malloc( sizeof(float)*iwl);





            for (int x=0 ; x < iWidth;  x++) {

                // reduce the size of the comparison window if we are near the boundary
                int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,MIN(iHeight-1-y,MIN(x,y))));


                // research zone depending on the boundary and the size of the window
                int imin=MAX(x-iDBloc,iDWin0);
                int jmin=MAX(y-iDBloc,iDWin0);

                int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);



                //  clear current denoised patch
                for (int ii=0; ii < iChannels; ii++) fpClear(fpODenoised[ii], 0.0f, iwl);



                // maximum of weights. Used for reference patch
                float fMaxWeight = 0.0f;


                // sum of weights
                float fTotalWeight = 0.0f;


                for (int j=jmin; j <= jmax; j++)
                    for (int i=imin ; i <= imax; i++)
                        if (i!=x || j!=y) {

                            float fDif = fiL2FloatDist2(fpI,fpI,x,y,i,j,iDWin0,iChannels,iWidth,iWidth);

                            // dif^2 - 2 * fSigma^2 * N      dif is not normalized
                            fDif = MAX(fDif - 2.0f * (float) icwl *  fSigma2, 0.0f);
                            fDif = fDif / fH2;

                            float fWeight = wxSLUT(fDif,fpLut);

                            if (fWeight > fMaxWeight) fMaxWeight = fWeight;

                            fTotalWeight += fWeight;


                            for (int is=-iDWin0; is <=iDWin0; is++) {
                                int aiindex = (iDWin+is) * ihwl + iDWin;
                                int ail = (j+is)*iWidth+i;

                                for (int ir=-iDWin0; ir <= iDWin0; ir++) {

                                    int iindex = aiindex + ir;
                                    int il= ail +ir;

                                    for (int ii=0; ii < iChannels; ii++)
                                        fpODenoised[ii][iindex] += fWeight * fpI[ii][il];

                                }

                            }


                        }




                // current patch with fMaxWeight
                for (int is=-iDWin0; is <=iDWin0; is++) {
                    int aiindex = (iDWin+is) * ihwl + iDWin;
                    int ail=(y+is)*iWidth+x;

                    for (int ir=-iDWin0; ir <= iDWin0; ir++) {

                        int iindex = aiindex + ir;
                        int il=ail+ir;

                        for (int ii=0; ii < iChannels; ii++)
                            fpODenoised[ii][iindex] += fMaxWeight * fpI[ii][il];

                    }
                }



                fTotalWeight += fMaxWeight;







                // normalize average value when fTotalweight is not near zero
                if (fTotalWeight > fTiny) {



                    for (int is=-iDWin0; is <=iDWin0; is++) {
                        int aiindex = (iDWin+is) * ihwl + iDWin;
                        int ail=(y+is)*iWidth+x;

                        for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                            int iindex = aiindex + ir;
                            int il=ail+ ir;

                            fpCount[il]++;

                            for (int ii=0; ii < iChannels; ii++) {
                                fpO[ii][il] += fpODenoised[ii][iindex] / fTotalWeight;

                            }

                        }
                    }


                }






            }



            for (int ii=0; ii < iChannels; ii++) free(fpODenoised[ii]);
            free(fpODenoised);


        }




    }









    for (int ii=0; ii < iwxh; ii++)
        if (fpCount[ii]>0.0) {
            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] /= fpCount[ii];

        }       else {

            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] = fpI[jj][ii];
        }




    // delete memory
    free(fpLut);
    free(fpCount);



}
