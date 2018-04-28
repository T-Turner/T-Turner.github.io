---
permalink: /ab_pde_heat
layout: default
---

### Homework Adams Bashforth PDE Heat Equation Solver
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves the heat equation with the initial boundary condition of U(x,0)=2 * x + sin(2 * pi*x) + 1. 

### Input:
The only input required is the boundary conditions, beta, total time, initial condtition, and length.

### Output: 
It will output the temperature at 9 locations at time is zero, and total time/10 and so o

### Usage:

```c++
  ab_pde();
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<iomanip>

//u_t=beta*u_xx
namespace
{
	const int n = 11, val = 100;//val needs to be divisible by 10 (bookkeeping purposes)
	double a2[n][n], u[n][n], b[n], k[5][n];
	double U[val][n], temp[n][n], temp1[n][n], temp2[n][n], utemp[n][n];

	double x[n], x1 = 0, x2 = 1, t = 0, time = 12000;
	double dx = (x2 - x1) / (n - 1);
	double dt = time / val;
	double beta = 0.00001;
	double s = beta * dt / (dx*dx);//gain
	const double pi = 3.1415926535;

	double low[n][n], up[n][n], z[n], y[n];
}

double f(double x)//IVP
{
	return 2 * x + sin(2 * pi*x) + 1;
}
void check_gain()
{
	if (s > 0.5)
		std::cout << "The number of steps specified leads to an unstable solution" << std::endl;
}
void boundsandx()
{
	for (int i = 0; i < n; ++i)
		x[i] = dx * i;
	b[0] = s;
	b[n - 1] = 4 * s*dx;
}
void spdiag()
{
	double diag = 1 + 2 * s;
	double subsup = -s;
	for (int j = 0; j < n; ++j)
	{
		a2[j][j] = diag;
		if (j + 1 < n)
			a2[j][j + 1] = subsup;
		if (j - 1 >= 0)
			a2[j][j - 1] = subsup;
	}
	a2[n - 1][n - 2] = -2 * s;
}

void matadd(double mat1[n][n], double mat2[n], double result[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
		{
			if (col == 0)
				result[row][col] = mat1[row][col] + mat2[row];
			else
				result[row][col] = mat1[row][col];
		}
}

void lu_tri(double A[n][n])
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
void find_y(double B[n][n])
{
	y[0] = B[0][0] / low[0][0];
	for (int i = 1; i < n; ++i)
	{
		for (int k = 0; k < i + 0; ++k)
			y[i] += low[i][k] * y[k];
		y[i] = (B[i][0] - y[i]) / low[i][i];
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
void lu_decom(double B[n][n])
{
	find_y(B);
	find_z();
}

void print(double mat[val][n])
{
	for (int i = 0; i < val; i += val / 10)
	{
		for (int j = 0; j < n; ++j)
			std::cout << mat[i][j] << " ";
		std::cout << std::endl;
	}
}

void cleanup()
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			temp[i][j] = 0;
			temp1[i][j] = 0;
			temp2[i][j] = 0;
			z[i] = 0;
			y[i] = 0;
		}
}

void initial()
{
	check_gain();//check stability
	boundsandx();//fills boundary, x arrays
	spdiag(); //creates the diagonalized matrix
	lu_tri(a2);
	for (int i = 0; i < n; ++i)
		u[0][i] = f(x[i]);
	for (int i = 0; i < n; ++i)
		U[0][i] = u[0][i];
}

void emulate(double orig[n][n], double emul[n])
{
	for (int i = 0; i < n; ++i)
		orig[0][i] = emul[n];
}
void setzero(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = 0;
}
void addk(double matu[n][n], int knum, int j)
{
	setzero(utemp);
	for (int i = 0; i < n; ++i)
		utemp[0][i] = U[j-1][i] + k[knum][i] * .5;
}

void logic(double matu[n][n], int knum, int j)
{
	addk(matu, knum, j);
	cleanup();
	for (int i = 0; i < n; ++i)
		temp[i][0] = utemp[0][i];
	matadd(temp, b, temp2);
	lu_decom(temp2);
	for (int i = 0; i < n; ++i)
		k[knum][i] = z[i];
}

void ab_pde()
{
	initial();
	for (int j = 1; j < val; ++j)
	{
		logic(u, 1, j);
		logic(utemp, 2, j);
		logic(utemp, 3, j);
		logic(utemp, 4, j);
		for (int i = 0; i < n; ++i)
			U[j][i] = k[1][i] / 6 + k[2][i] / 3 + k[3][i] / 3 + k[4][i] / 6;
		for (int i = 1; i < 5; ++i)
			for (int j = 0; j < n; ++j)
				k[i][j] = 0;
	}
	print(U);
}

int main()
{
	ab_pde();
	getchar();
	return EXIT_SUCCESS;
}

```


### Last Modified:
April 28, 2018
