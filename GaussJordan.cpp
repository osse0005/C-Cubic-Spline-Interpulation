/**************************************************************************
 Filename:			GaussJordan.cpp
 Version:			1.0
 Author:			Andrew Tyler MODIFIED
 Student #:			040 585 009
 Course Name:		Numerical Analysis CST8233
 Lab Sect:			310
 Assignment #:		3
 Assignment Name:	Cubic Spline
 Due Date:			December 8th, 2012
 Submission Date:	December 8th 2012
 Professor:			Andrew Tyler
 Purpose:			Preform Gauss Jordan algorthim
 ***************************************************************************/

#include "GaussJordan.h"
/**************************************************************************
 Function name:
 Purpose:
 In parameter:
 Out parameter:
 Version:		1
 Author:		Andrew Tyler (Notes)
 ***************************************************************************/
void gaussj(float** a, int n, float** b, int m){
    int i, icol, irow, j, k, l, ll;
    float pivinv, big, dum, temp;
    
    /*temporary storage on the heap to record which rows have been reproduced (to the unit matrix)*/
    int* ipiv=(int*)malloc(n*sizeof(int));
    
    /*initialize elements to 0*/
    for(j=0;j<n;j++) ipiv[j]=0;
    
    /*the outer loop for all the rows to reduce*/
    for(i=0;i<n;i++){
        /*initially we have no plot*/
        big=0.0;
        /*search all the rows*/
        for(j=0;j<n;j++){
            /*but not those that have already been reduced*/
            if(ipiv[j] != 1){
                /*search all the columns of this row*/
                for(k=0;k<n;k++){
                    /*expect 0 if this row has not been previously reduced*/
                    if(ipiv[k] == 0){
                        /*record this pivot if it is the current max (j=row, k=col*/
                        if (fabs(a[j][k]) >= big){
                            big = fabs (a[j][k]);
                            /*the row its in*/
                            irow = j;
                            /*the columns its in*/
                            icol = k;
                        }/*end if (fabs)*/
                        /*a column can only be reduced once - if more than that, something's broken*/
                    }/*end if (k)*/else if(ipiv[k] > 1)printf("gaussj: Singular Matrix-1");
                }/*end for (k)*/
            }/*end if (j)*/
        }/*end for (j)*/
        
        /*record that we found a pivot for row icol*/
        ++(ipiv[icol]);
        
        /*if irow!=icol, need to swop rows to bring the pivot to the diagonal*/
        if (irow != icol){
            for(l=0;l<n;l++) SWAP(a[irow][l], a[icol][l]);
            for(l=0;l<m;l++) SWAP(b[irow][l], b[icol][l]);
        }/*end if*/
        
        
        /*if the pivot is 0 then we'll be dividing by 0*/
        if (a[icol][icol] == 0.0) printf("gaussj: Singular Matrix-2");
        
        /*convert the diagonal to 1 and adjust the row in proportion*/
        pivinv=1.0/a[icol][icol];
        for(l=0;l<n;l++) a[icol][l] *= pivinv;
        for(l=0;l<m;l++) b[icol][l] *= pivinv;
        
        /*convert rest of column to 0 and adjust the rows accordingly*/
        for(ll=0;ll<n;ll++)
            if(ll != icol){
                dum=a[ll][icol];
                for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
            }/*end if*/
    }
    
    free(ipiv);
}/*end gaussj*/
