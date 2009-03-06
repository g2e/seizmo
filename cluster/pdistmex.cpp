/*
 * pdistmex.cpp
 *
 * Calculates pairwise distances between observations.
 * Helper function to pdist.m
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1993-2006 The MathWorks, Inc.
 */

/* $Revision: 1.1.6.5 $  $Date: 2007/03/15 19:27:03 $ */

#include "mex.h"
#include <math.h>
#include <string.h>

/* Euclidean distance */
template<class T>
void eucdist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
     */
    mwIndex i,j,k;
    T   theSum,Y;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y;
            }
            *(d++) = (T)sqrt(theSum);
        }
    }
}

/* Standardized Euclidean distance */
template<class T>
void seudist(T *x, mwSize m, mwSize n, T *arg, T *d)
{
    /*
     d = sqrt(((XI-XJ).^2) * arg);           % Standardized Euclidean
     */
    /* arg is a column vector of 1/var(X) */
    mwIndex i,j,k;
    T   theSum,Y;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += Y*Y*(*(arg+k));
            }
            *(d++) = (T)sqrt(theSum);
        }
    }
}

/* City Block Distance */
template<class T>
void citdist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     d = sum(abs((XI-XJ)),2);                % City Block
     */
    mwIndex i,j,k;
    T   theSum,Y;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (*XI)-(*XJ);
                theSum += (T)fabs(Y);
            }
            *(d++) = theSum;
        }
    }
}

/* Mahalanobis distance */
template<class T>
void mahdist(T *x, mwSize m, mwSize n, T *arg, T *d)
{
    /*
     Y = XI - XJ;
     d = sqrt(sum((Y*arg).*Y,2));            % Mahalanobis
     */
    /* arg is inv(cov(X)) */
    mwIndex i,j,k,l;
    T   theSum, inner,Y,YY;
    T   *XXI, *XXJ;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;

            for (l=0;l<n; l++,XI++,XJ++){
                XXJ = x+j*n;
                XXI = x+i*n;
                inner = 0;
                for (k=0;k<n;k++,XXI++,XXJ++){
                    YY = (*XXI)-(*XXJ);
                    inner += YY*arg[k+l*n];
                }
                Y = (*XI)-(*XJ);
                theSum += inner * Y;
            }
            *(d++) = (T)sqrt(theSum);
        }
    }
}

/* Minkowski distance */
template<class T>
void mindist(T *x, mwSize m, mwSize n, T arg, T *d)
{
    /*
     d = sum(abs((XI-XJ)).^arg,2).^(1/arg);  % Minkowski
     */
    mwIndex i,j,k;
    T   theSum,Y;
    T   argRecip = 1/arg;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                Y = (T)fabs((*XI)-(*XJ));
                theSum += (T)pow(Y,arg);
            }
            *(d++) = (T)pow(theSum,argRecip);
        }
    }
}

/* Cosine and Correlation distances */
template<class T>
void coscordist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     d = 1 - sum(XI.*XJ,2);   % Cosine & Corr & RankCorr
     */

    /* This actually calculates the dot product of pairs of vectors.  It
     * assumes that the data have been properly preprocessed:  ranked for
     * Spearman's rank correlation distance, normalized to zero mean for
     * both linear and rank correlation distances, and normalized to unit
     * length for all three distances.
     */
    mwIndex i,j,k;
    T   theSum;
    T   *XI, *XJ, *XI0;

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            theSum = 0;
            for (k=0;k<n;k++,XI++,XJ++){
                theSum += (*XI)*(*XJ);
            }
            /* theSum may overshoot 1 due to round-off, protect against
             * that without overwriting NaNs */
            *(d++) = theSum>1 ? 0 : 1-theSum;
        }
    }
}

/* Hamming distance */
template<class T>
void hamdist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     d = sum(XI ~= XJ,2) / size(XI,2);       % Hamming
     */
    mwIndex i,j,k;
    T   theSum;
    T   *XI, *XJ, *XI0;
    T   *theNaN = (T*)mxCalloc(m,sizeof(T));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI)!=(*XJ)) {
                        theSum++;
                    }
                }
                *(d++) = theSum/n;
            }
        }
    }
    mxFree(theNaN);
}

/* Jaccard distance */
template<class T>
void jacdist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     nz = XI ~= 0 | XJ ~= 0;
     ne = XI ~= XJ;
     d = sum(ne&nz,2) ./ sum(nz,2);          % Jaccard
     */
    mwIndex i,j,k;
    T   theSum,nz;
    T   *XI, *XJ, *XI0;
    T   *theNaN = (T*)mxCalloc(m,sizeof(T));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theSum = 0;
                nz = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    if ((*XI) || (*XJ)) {
                        nz++;
                        if ((*XI)!=(*XJ)) {
                            theSum++;
                        }
                    }
                }
                if (nz) {
                    *(d++) = theSum/nz;
                } else {
                    *(d++) = (T)mxGetNaN();
                }
            }
        }
    }
    mxFree(theNaN);
}

/* Chebychev distance */
template<class T>
void chedist(T *x, mwSize m, mwSize n, T *d)
{
    /*
     d = max(abs(XI-XJ),[],2);
     */
    mwIndex i,j,k;
    T   theMax,Y;
    T   *XI, *XJ, *XI0;
    T   *theNaN = (T*)mxCalloc(m,sizeof(T));

    /* Determine which rows of X have missing data.
     */
    XI = x;
    for (i=0; i<m; i++) {
        for (k=0;k<n;k++,XI++) if (mxIsNaN(*XI)) theNaN[i] = *XI;
    }

    XI = x;
    for (i=0; i<m; i++) {
        XI0 = XI;
        XJ = XI+n;
        for (j=i+1; j<m; j++) {
        /* XI = x + i*n; XJ = x + j*n; */
            XI = XI0;
            /* If XI or XJ have missing data, set their distance to NaN.
             */
            if (theNaN[i] || theNaN[j]) {
                XI += n;
                XJ += n;
                *(d++) = theNaN[i] + theNaN[j];
            } else {
                theMax = 0;
                for (k=0;k<n;k++,XI++,XJ++){
                    Y = (T)fabs((*XI)-(*XJ));
                    if (Y>theMax) {
                        theMax = Y;
                    }
                }
                *(d++) = theMax;
            }
        }
    }
    mxFree(theNaN);
}


/* the dispatcher function */
template<class T>
void distfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], T classDummy)
{
    int     status;
    mwSize  numCoords,numPoints;
    char    metric[4];
    T       *x,*d,*arg,scalarArg;

    /*  get the metric */
    status = mxGetString(prhs[1],metric,4);

  /*  create a pointer to the input matrix y */
    x = (T*)mxGetData(prhs[0]);

  /*  get the dimensions of the matrix input y */
    numCoords = mxGetM(prhs[0]);
    numPoints = mxGetN(prhs[0]);

  /* get extra arg  */
    if (nrhs>2 && !mxIsEmpty(prhs[2])) {
        if (strcmp(metric,"min") == 0) {  /*minkowski takes a scalar of any type */
            scalarArg = (T)mxGetScalar(prhs[2]);
        } else { /*seucidean takes vector, mahalanobis takes matrix, type must match */
            if (mxGetClassID(prhs[2]) == mxGetClassID(prhs[0])) {
                arg = (T*)mxGetData(prhs[2]);
            } else {
                mexErrMsgIdAndTxt("stats:pdistmex:MixedInputTypes",
                                  "Additional input arguments must be the same class as X.");
            }
        }
    }

    /* make sure that the distance matrix can be created, then create it.  doing
     * this in double remains exact except in the cases where we error out anyway. */
    double numDists = ((double)numPoints * (double)(numPoints-1)) / 2;
    if (numDists >= (double)MWSIZE_MAX) {
        mexErrMsgIdAndTxt("stats:pdistmex:OutputTooLarge",
                          "Distance matrix has more elements than the maximum allowed size in MATLAB.");
    }
    plhs[0] = mxCreateNumericMatrix(1, (mwSize)numDists, mxGetClassID(prhs[0]), mxREAL);
    
  /*  create a pointer to a copy of the output matrix */
    d = (T*)mxGetData(plhs[0]);

  /*  call the appropriate distance subroutine */
    if (strcmp(metric,"euc") == 0)
        eucdist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"seu") == 0)
        seudist(x,numPoints,numCoords,arg,d);
    else if(strcmp(metric,"cit") == 0)
        citdist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"min") == 0)
        mindist(x,numPoints,numCoords,scalarArg,d);
    else if(strcmp(metric,"cos") == 0)
        coscordist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"cor") == 0)
        coscordist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"spe") == 0)
        coscordist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"ham") == 0)
        hamdist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"jac") == 0)
        jacdist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"che") == 0)
        chedist(x,numPoints,numCoords,d);
    else if(strcmp(metric,"mah") == 0)
        mahdist(x,numPoints,numCoords,arg,d);
}

/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /*  check for proper number of arguments */
    if (nrhs<2) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooFewInputs",
                          "Two input arguments required.");
    } else if(nlhs>1) {
        mexErrMsgIdAndTxt("stats:pdistmex:TooManyOutputs",
                          "Too many output arguments.");
    }

  /* Check the type of the input array */
  /* Currently only works with double or single(float) */
    if (mxIsDouble(prhs[0])) {
        distfun(nlhs, plhs, nrhs, prhs, (double)1.0);
    } else if (mxIsSingle(prhs[0])) {
        distfun(nlhs, plhs, nrhs, prhs, (float)1.0);
    } else {
        mexErrMsgIdAndTxt("stats:pdistmex:BadInputType",
                          "PDISTMEX only supports real DOUBLE and SINGLE data.");
    }
}
