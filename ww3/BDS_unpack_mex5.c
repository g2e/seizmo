/*

BDS_unpack_mex5.c
GRiB Binary Data Section decoding
Relevant code from W. Ebisuzaki, NCEP

30 Aug, 2005
   fixed alloc sizeof bug noted by several users under MATLAB7/R14
   Both Felipe Nievinski and Julien Choisnard noticed that seg viols
   occured because of the sizeof(mxREAL).  These have been changed to
   size(float).  

05 Feb 2010
   replaced UNDEFINED=10E20 with NaN. added a few comments
*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

/* PROTOTYPES */
void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
        int n_bits, int n, double ref, double scale, double );
	
/* undefined value -- if bitmap */
#define UNDEFINED		9.999e20 

/************************************************************

  ####     ##     #####  ######  #    #    ##     #   #
 #    #   #  #      #    #       #    #   #  #     # #
 #       #    #     #    #####   #    #  #    #     #
 #  ###  ######     #    #       # ## #  ######     #
 #    #  #    #     #    #       ##  ##  #    #     #
  ####   #    #     #    ######  #    #  #    #     #

************************************************************/

void mexFunction(int            nlhs,
                 mxArray      *plhs[],
                 int            nrhs,
                 const mxArray *prhs[])
{
   unsigned char *databits, *bitmap;

   int n_bits,n,i,mb,nb;
   double ref,scale,*tempd;
   FILE *tempfile,*fopen();
   float *tempf;
   double NaN=mxGetNaN();

   bool bms_empty;
   
   databits=(unsigned char *)mxGetData(prhs[0]);   /* bds_struct.bindata  */
   bms_empty=mxIsEmpty(prhs[1]);                   /* check for bitmap in bms_struct.bitmap */

   if(bms_empty)
      bitmap=NULL;
   else
      bitmap=(unsigned char *)mxGetData(prhs[1]);  /* get bitmap, if it's there    */
   
   n_bits=mxGetScalar(prhs[2]);                    /* bds_struct.nbits             */
   n=mxGetScalar(prhs[3]);                         /* nxny                         */
   ref=mxGetScalar(prhs[4]);                       /* decimal_sf*bds_struct.RefVal */
   scale=mxGetScalar(prhs[5]);                     /* decimal_sf*bin_sf            */
   
   /* tempf=(float *)mxMalloc(sizeof(mxREAL)*n); */
   /* sizeof(mxREAL) caused seg viol in R14 */
   tempf=(float *)mxMalloc(sizeof(float)*n);

   /* mb=mxGetM(prhs[0]); */
   /* nb=mxGetN(prhs[0]); */

   /* call unpack routine 
      tempf is returned; it is the decoded binary data 
      NaN is used as the undef value */
   BDS_unpack(tempf, databits, bitmap, n_bits, n, ref, scale, NaN);
   
   plhs[0]=mxCreateDoubleMatrix(n,1,mxREAL);
   tempd=mxGetPr(plhs[0]);
   for (i=0;i<n;i++){ tempd[i]=tempf[i]; } 
   mxSetPr(plhs[0],tempd);
   
   return;
}

/* for simple unpacking of a grid */
/* wesley ebisuzaki, NCEP, wgrib */
/* http://wesley.wwb.noaa.gov */

void BDS_unpack(float *flt, 
                unsigned char *bits,        /* bds_struct.bindata  */  
		unsigned char *bitmap,      /* bms_struct.bitmap */
                int n_bits,                 /* bds_struct.nbits             */
                int n,                      /* nxny                         */
                double ref,                 /* decimal_sf*bds_struct.RefVal */
                double scale,               /* decimal_sf*bin_sf            */
		double FillVal)
		{

    int i, j, k;
    unsigned int map_mask, bit_mask;

    map_mask = bit_mask = 128;

    for (i = 0; i < n; i++) {
	if (bitmap) {
	    j = (*bitmap & map_mask);
	    if ((map_mask >>= 1) == 0) {
		map_mask = 128;
		bitmap++;
	    }
	    if (j == 0) {
		*flt++ = FillVal;
		continue;
	    }
	}

	j = 0;
	k = n_bits;
	while (k) {
	    if (k >= 8 && bit_mask == 128) {
		j = 256 * j + *bits;
		bits++;
		k -= 8;
	    }
	    else {
	        j = j + j + ((*bits & bit_mask) != 0);
		if ((bit_mask >>= 1) == 0) {
		    bits++;
		    bit_mask = 128;
		}
		k--;
	    }
	}
	*flt++ = ref + scale*j;
   }
   return;
}
