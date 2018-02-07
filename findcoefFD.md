---
permalink: /findcoefFD
layout: default
---

### Name: findcoefFD
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program finds the coefficients of a finite difference approximation using a LU decomposition.

### Input:
The program requires a value for k, the real value, and the estimated or calulated value.

### Output: 
The progam will output the coefficients required for finite difference approximations

### Usage:
This program requires that it be called in the following fashion:
```c++
  findcoef(k, xbar, x)
```
### Implementation/Code:
The code is implemented in the following manner:
```c++
#include<iostream>
#include<iomanip>

namespace
{
	const int n;//must initialize n  here.
	double A[n][n], low[n][n], up[n][n], B[n], z[n], y[n];
}

void print(double mat[])
{
	std::cout << std::setprecision(3);
	for (int i = 0; i < n; ++i)
		std::cout << mat[i] << std::endl;
}

void lu_tri()
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			low[i][0] = A[i][0];
			up[i][j] = A[i][j];
		}
	for (int k = 0; k < n; ++k)
	{
		for (int col = n - 1; col >= 0; --col)
			up[k][col] = up[k][col] / up[k][k];
		for (int row = k + 1; row < n; ++row)
		{
			for (int col = n - 1; col >= k; --col)
				up[row][col] = up[row][col] - up[row][k] * up[k][col];
			low[row][k + 1] = up[row][k + 1];
		}
	}
}

void find_y()
{
	y[0] = B[0] / low[0][0];
	for (int i = 1; i < n; ++i)
	{
		for (int k = 0; k < i + 0; ++k)
		{
			y[i] += low[i][k] * y[k];
			//std::cout << " y " << y[i];
		}
		y[i] = (B[i] - y[i]) / low[i][i];
	}
}

void find_z()
{
	z[n - 1] = y[n - 1] / up[n - 1][n - 1];
	for (int i = n - 2; i >= 0; --i)
	{
		for (int k = i + 1; k < n; ++k)
			z[i] -= up[i][k] * z[k];
		z[i] = (y[i] + z[i]) / up[i][i];
	}
}

void lu_decom()
{
	lu_tri();
	find_y();
	find_z();
}

void chngmat(double mat[n][n], int num)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = num;
}

int factorial(int num)
{
	if (num > 1)
		return num * factorial(num - 1);
	else
		return 1;
}

void matprint(double mat[n][n])
{
	std::cout << std::setprecision(3);
	for (int row = 0; row < n; ++row)
	{
		for (int col = 0; col < n; ++col)
			std::cout << mat[row][col] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void findcoef(int k, double xbar, double x)
{
	chngmat(A, 1);
	for (int i = 1; i < n; ++i)
		for (int j = 0; j < n; ++j)
			A[i][j] = pow((x - xbar), i) / factorial(i - 1);
	for (int i = 0; i < n; ++i)
		B[i] = 0;
	B[k + 1] = 1;
	lu_decom();
	print(z);
}
```
### Output
A array of the coefficients in the z array.

### Last Modified:
February 7, 2018
