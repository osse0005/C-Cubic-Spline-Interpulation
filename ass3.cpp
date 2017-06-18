/**************************************************************************
 Filename:			ass3.cpp
 Version:			1.1.1
 Author:			Tim Osse
 Student #:			040 585 009
 Course Name:		Numerical Analysis CST8233
 Lab Sect:			310
 Assignment #:		3
 Assignment Name:	Cubic Spline
 Due Date:			December 8th, 2012
 Submission Date:	December 8th 2012
 Professor:			Andrew Tyler
 Purpose:			Calculate cubic spline interpulation, and exponential
					trend using least squares linear regression 
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include "GaussJordan.h"

using namespace std;

char* fileName;//file name to open
FILE* fp;//file pointer
int numPts = 0;
float b, a, c2;
float b_b;
float a_a;
float* x;//array of x values
float* y;//array of y values
float* yy;//array of y values


/********************************************************************
Function name:	ReadFile
Purpose:		read in x and y values from file into memory
In parameters:	none
Out parameters: void
Version:		1.2.3
Author:		    Tim Osse
**********************************************************************/
void ReadFile(void){

	float temp_x = 0.0f;
	float temp_y = 0.0f;
	int i,j = 0;
	int k;

	printf("OPEN A FILE OF DATA TO SPLINE...");
	fileName = (char*)malloc(sizeof(char)*25);
	printf("\nPlease enter the file name: ");
	scanf("%s",fileName);

	fp = fopen(fileName,"r");

	if(fp == NULL){
		printf("Empty file\n\n");
	}
	else {
		numPts=-1;
		while((k = fgetc(fp)) != EOF){
			if(k == '\n'){
				numPts++;
			}
		}
		//allocate memory for arrays
		x = (float*) malloc(sizeof(float)*numPts);
		y = (float*) malloc(sizeof(float)*numPts);
		yy = (float*) malloc(sizeof(float)*numPts);

		printf("number of data points = %d\n",numPts);
		rewind(fp);
		while(fscanf(fp,"%f %f",&temp_x,&temp_y) != EOF){
			x[j] = temp_x;
			y[j] = temp_y;
			yy[j] = temp_y;
			j++;
		}

	}//end else

	free(fileName);
}//end func


/**************************************************************************
 Function name: SplineTridiagonal
 Purpose:
 In parameter: int n
 Out parameter: A
 Version:		1
 Author:		Andrew Tyler (Notes)
 ***************************************************************************/
float** SplineTridiagonal (int n){
    
    int i, j, k = -1;
    float pattern [] = {1,4,1};
    
    float** A = new float *[n];
    for (i=0; i<n; i++) A[i] = new float[n];
    for (i=0;i<n;i++)      /*initialize all elements to 0*/
        for (j=0;j<n;j++)
            A[i][j]=0;
    
    /*put in the numbers*/
    A[0][0] = A[n-1] [n-1] = 1;     /*the first and last rows*/
    for (i=1;i<(n-1);i++)           /*all the rest*/
        memcpy((A[i]+(++k)), pattern, 3*sizeof(float));
    
    return A;
    
}

/**************************************************************************
 Function name: Splineb
 Purpose:
 In parameter: float* x, float* y, int n
 Out parameter: b
 Version:		1
 Author:		Andrew Tyler (Notes)
 ***************************************************************************/
float** Splineb(float* x, float* y, int n){
    
    int i, j;
    float h = fabs(x[1]-x[0]);        /*assume uniformly spaced data*/
    
    /*allocate space on the heap for n pointers*/
    
    float** b = new float*[n];  /*1 column but it's still a float** */
    
    /*each points to a single float*/
    for (i=0;i<n;i++) b[i] = new float[1];
    for (i=0;i<n;i++)       /*init all elements to 0*/
        for(j=0;j<1;j++)
            b[i][j]=0;
    
    /*calculate and set the elements of the matrix*/
    b[0][0] = b[n-1][0] = 0;    /*the first and last are 0*/
    
    for (i=1;i<(n-1);i++){
        b[i][0] = (y[i-1]-2*y[i]+y[i+1])*(6/(h*h));
    }/*end for*/
    
    return b;
    
}

/**************************************************************************
 Function name: splineAT
 Purpose:
 In parameter: float ax[], float ya[], float** y2a, int n, float x
 Out parameter: a*dx*dx*dx+b*dx*dx+c*dx+d;
 Version:		1
 Author:		Andrew Tyler (Notes)
 ***************************************************************************/
float splintAT(float xa[], float ya[], float** y2a, int n, float x){
    int klo = 0, khi, k;
    float h;
    
    /*Find the right place in the table by means of bisection*/
    khi = n-1;
    while (khi-klo > 1) {
        k = (khi+klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo = k;
    }/*end while*/
    
    /*klo and khi now bracket the input value of x*/
    h=xa[khi]-xa[klo];
    if(h==0.0) printf("Bad xa input to routine splint");
    
    /*do the interpolation*/
    float a,b,c,d;
    
    a=(y2a[khi][0]-y2a[klo][0])/(6*h);
    b=y2a[klo][0]/2;
    c=((ya[khi]-ya[klo])/h) - ((y2a[khi][0]+2*y2a[klo][0])*h/6);
    d=ya[klo];
    
    float dx=x-xa[klo];
    return a*dx*dx*dx+b*dx*dx+c*dx+d;
    
}

/********************************************************************
Function name:	LinearRegression
Purpose:		Calculate least squares fitting of data 
In parameters:	float* x, float* y, int ndata, float* m, float* c, float* sigm, float*
sigc, float* sigma
Out parameters: *a
Version:		1.1.2
Author:		    Tim Osse
**********************************************************************/
void linearRegression(float* x, float* y, int ndata, float* b, float* a)
{
	int i;
	float t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat,sumSquares=0.0,sigma2,residual;
	*b=0.0;

	for(i=0;i<ndata;i++){
		y[i] = log(y[i]);
	}

	for(i=0;i<ndata;i++){
		sx += x[i];
		sy += y[i];
	}
	ss=ndata;
	sxoss=sx/ss;
	for(i=0;i<ndata;i++){
		t=x[i]-sxoss;
		st2 += t*t;
		*b += t*y[i];
	}
	*b /= st2;				//Solve for b
	*a=(sy-sx*(*b))/ss;		//a


	*a = exp(*a);
}

/**************************************************************************
 Function name: main
 Purpose: main menu, call functions, print results
 In parameter: n/a
 Out parameter: n/a
 Version:		1.1.2
 Author:		Tim Osse
 ***************************************************************************/
int main(void){

	float** A;
	float** b;
	float interpolate;
	float exponential;
	int response;
	float i;
	float increment;
	float range;

	while(1){
		cout<<"\nMENU";
		cout<<"\n1: Read in a new file";
		cout<<"\n0: Quit:\n";
		cin>>response;
		if(response == 0)return 0;

		ReadFile();
		A = SplineTridiagonal(numPts);//make the tridiagonal matrix
		b=Splineb(x,yy,numPts);
		gaussj(A, numPts, b, 1);
		linearRegression(x,y,numPts,&b_b,&a);

		cout<<setprecision(3)<<"Exponential Fit is: y = ("<<a<<"*exp("<<b_b<<"*x)";
		cout<<"\nPlease enter the range to interpolate to (cannot exceed 5): ";
		cin>>range;
		cout<<" "<<endl;
		if(range > 5){cout<<"INVALID INPUT"<<endl; return 0;}//check range value
		increment = range / 10;

		for(i =0;i<=range;i+=increment){

			interpolate = splintAT(x, yy, b, numPts, i);
			exponential = a*exp(b_b*i);

			std::cout << std::fixed << std::setprecision(2);
			cout<<"x = "<<i<<setw(20)<<"interpulation = "<<interpolate<<setw(20)<<"exponential = "<<exponential<<endl<<endl;		//publish the interpulated y

		}
	}
	free(x);
	free(y);
	free(yy);
	return 0;
}

