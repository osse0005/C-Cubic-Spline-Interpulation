/**************************************************************************
 Filename:			GaussJordan.h
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
 Purpose:			header for Gauss Jordan class
 ***************************************************************************/

#ifndef _GAUSSJORDAN_H_
#define _GAUSSJORDAN_H_

#include <math.h>
#include <iostream>

#define SWAP(a,b) {temp =(a);(a) = (b);(b)=temp;}

void gaussj(float** , int , float** , int );

#endif //_GAUSSJORDAN_H_