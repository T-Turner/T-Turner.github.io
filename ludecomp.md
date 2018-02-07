---
permalink: /ludecomp
layout: default
---

## Initial Conditions

### Name: initialcond
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program used LU decomposition to solve a system of linear equations. The relationship of nxn matrix to time is dramatic, its path log curve, but can be approximated using a 2nd order polynomial time=5E-7xn +0.0001*n.

### Input:
The program requires the matrix A, and B be defined in the namespace.
 
### Output: 
This program outputs the solution in the z array.

### Usage:
```c++
  lu_decom();
```
### Implementation/Code:
These programs are implemented in the following manner. 
```c++
#include<iostream>
#include<iomanip>

namespace
{
	const int n = 5;
	double A[n][n];//insert A matrix = {row1.1,row1.2...,row1.n,row2.1-row2.n---rown.n};
	double B[n];//insert B matrix ={a,b,.....}
	double low[n][n], up[n][n], z[n], y[n];
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
			y[i] += low[i][k] * y[k];
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
	print(z);
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

int main()
{
	matprint(A);
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
February 7, 2018
