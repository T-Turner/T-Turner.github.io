---
permalink: /5p_laplacian
layout: default
---

## Five Point Laplacian
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
	const int n = 5; const int m = n - 1;
	const int N = m * m - 2 * m + 2;
	const int stuff = N - 1;
	double fill[n][n], a[N][N], b[N];
	double A[stuff][stuff], B[stuff];
	double top = 6, bot = 8, rside = 5, lside = 7;
	double low[stuff][stuff], up[stuff][stuff], z[stuff], y[stuff];
}

void fill_bound(double mat[n][n])//no corners
{
	for (int i = 1; i < m; ++i)
	{
		mat[0][i] = top;
		mat[m][i] = bot;
		mat[i][0] = lside;
		mat[i][m] = rside;
		for (int j = 1; j < m; ++j)
			mat[i][j] = 1;
	}
}

void five_ps(double mata[N][N], double matb[N], double matf[n][n])
{
	int k = 0;
	for (int r = 1; r < m; ++r)
		for (int c = 1; c < m; ++c)
		{
			if (r == 1)
				k = r + c - 1;
			if (r > 1)
				k = c + (r - 1)*(m - 1);

			if ((r - 1) == 0)
				matb[k] -= matf[r - 1][c];
			if ((r - 1) != 0)
				mata[k][k - m + 1] = -matf[r - 1][c];

			if ((r + 1) == m)
				matb[k] -= matf[r + 1][c];
			if ((r + 1) != m)
				mata[k][k + m - 1] = -matf[r + 1][c];

			if ((c - 1) == 0)
				matb[k] -= matf[r][c - 1];
			if ((c - 1) != 0)
				mata[k][k - 1] = -matf[r][c - 1];

			if ((c + 1) == m)
				matb[k] -= matf[r][c + 1];
			if ((c + 1) != m)
				mata[k][k + 1] = -matf[r][c + 1];

			mata[k][k] = 4 * matf[r][c];
		}
}

void clean(double mata[N][N], double matb[N], double newa[stuff][stuff], double newb[stuff])
{
	for (int i = 1; i < N; ++i)
	{
		newb[i - 1] = matb[i];
		for (int j = 1; j < N; ++j)
			newa[i - 1][j - 1] = mata[i][j];
	}
}

void setup(double mata[N][N], double matb[N], double matf[n][n])
{
	fill_bound(matf);
	five_ps(mata, matb, matf);
}

void print(double mat[])
{
	for (int i = 0; i < stuff; ++i)
		std::cout << mat[i] << std::endl;
}

void lu_tri(double mata[stuff][stuff])
{
	for (int i = 0; i < stuff; ++i)
		for (int j = 0; j < stuff; ++j)
		{
			low[i][0] = mata[i][0];
			up[i][j] = mata[i][j];
		}
	for (int k = 0; k < stuff; ++k)
	{
		for (int col = stuff - 1; col >= 0; --col)
			up[k][col] = up[k][col] / up[k][k];
		for (int row = k + 1; row < stuff; ++row)
		{
			for (int col = stuff - 1; col >= k; --col)
				up[row][col] = up[row][col] - up[row][k] * up[k][col];
			low[row][k + 1] = up[row][k + 1];
		}
	}
}

void find_y(double matb[stuff])
{
	y[0] = matb[0] / low[0][0];
	for (int i = 1; i < stuff; ++i)
	{
		for (int k = 0; k < i + 0; ++k)
			y[i] += low[i][k] * y[k];
		y[i] = (matb[i] - y[i]) / low[i][i];
	}
}

void find_z()
{
	z[stuff - 1] = y[stuff - 1] / up[stuff - 1][stuff - 1];
	for (int i = stuff - 2; i >= 0; --i)
	{
		for (int k = i + 1; k < stuff; ++k)
			z[i] -= up[i][k] * z[k];
		z[i] = (y[i] + z[i]) / up[i][i];
	}
}

void lu_decom(double newa[stuff][stuff], double newb[stuff])
{
	lu_tri(newa);
	find_y(newb);
	find_z();
	print(z);
}

int main()
{
	setup(a, b, fill);
	clean(a, b, A, B);
	lu_decom(A, B);
	
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified
March 20, 2018
