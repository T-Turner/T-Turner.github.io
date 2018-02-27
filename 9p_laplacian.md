---
permalink: /9p_laplacian
layout: default
---

## Nine Point Laplacian
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program computes a mesh for a 5 point mesh and then solves the linear system.

### Input:
function for the laplacian. once the function is entered, the boundary conditions are needed.

### Output: 
The solution to the ODE. 

### Usage:
This program requires that it be called in the following fashion:
```c++
five_stencil(a);//a is the nxn matrix
```

### Implementation/Code:
The code is implemented in the following manner:
```c++
#include<iostream>

namespace
{
	const int n = 9, m = n - 1;
	double a[n][n];//boundaries are row and columns 0, & n-1
	double low[n][n], up[n][n], z[n], y[n], b[n][n];
}

void five_stencil(double mata[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			mata[i][i] = -4.0;
			if (i != m) {
				mata[i][i + 1] = 1.0;
				if ((i + 1) % 3 == 0)
					mata[i][i + 1] = 0.0;
			}
			if (i != 0) {
				mata[i][i - 1] = 1.0;
				if ((i + 1) % 3 == 1)
					mata[i][i - 1] = 0.0;
			}
			if (j == m - 3 && i <= m - 3)
				mata[i][i + 3] = 1.0;
			if (j == 3 && i >= 3)
				mata[i][i - 3] = 1.0;
		}
}

void five_b(double matb[n][n])
{
	for (int i = 0; i < n; ++i)
	{
		//b[i][0]=
	}
}

void lu_tri(double mata[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			low[i][0] = mata[i][0];
			up[i][j] = mata[i][j];
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

void find_y(double matb[n][n])
{
	y[0] = matb[0][0] / low[0][0];
	for (int i = 1; i < n; ++i)
	{
		for (int k = 0; k < i + 0; ++k)
			y[i] += low[i][k] * y[k];
		y[i] = (matb[i][0] - y[i]) / low[i][i];
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

void setequal(double matb[n][n])
{
	for (int i = 0; i < n; ++i)
		matb[i][0] = z[i];
}

void lu_decom(double mata[n][n], double matb[n][n])
{
	lu_tri(mata);
	find_y(matb);
	find_z();
	setequal(matb);
}

void print(double mat[n][n])
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			std::cout << mat[i][j] << " ";
		std::cout << std::endl;
	}
}

int main()
{
	five_stencil(a);
	print(a);
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified
February 26, 2018
