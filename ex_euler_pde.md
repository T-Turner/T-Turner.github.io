---
permalink: /ex_euler_pde
layout: default
---

### Explicit Euler for Solving Heat Equation
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves the heat equation

### Input:
The only input required is the boundary conditions, beta, total time, initial condtition, and length.

### Output: 
It will output the temperature at 9 locations at time is zero, and total time/10 and so on.

### Usage:

```c++
  euler();
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<iomanip>

//u_t=beta*u_xx
namespace
{
	const int n = 9, val = 600;//val needs to be divisible by 10 (bookkeeping purposes)
	double a2[n][n], u[n][n], b[n];
	double U[val][n];

	double x[n], x1 = 0, x2 = 1, t = 0, time = 12000;
	double dx = (x2 - x1) / (n + 1);
	double dt = time / val;
	double beta = 0.00001;
	double s = beta * dt / (dx*dx);//gain
	const double pi = 3.1415926535;
}
//u(x,0)=2x+sin(2*pi*x)+1
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
	b[n - 1] = 2 * s*dx;
}

void spdiag()
{
	double diag = 1 - 2 * s;
	double subsup = s;
	for (int j = 0; j < n; ++j)
	{
		a2[j][j] = diag;
		if (j + 1 < n)
			a2[j][j + 1] = subsup;
		if (j - 1 >= 0)
			a2[j][j - 1] = subsup;
	}
	a2[n - 1][n - 1] = 1 - s;
}

void matmult(double mat1[n][n], double mat2[n][n], double result[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			for (int last = 0; last < n; ++last)
				result[row][col] += mat1[row][last] * mat2[last][col];
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

void setzero(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = 0;
}
void setequal(double mat1[n][n], double mat2[n][n], int row)
{
	for (int i = 0; i < n; ++i)
		mat1[i][0] = mat2[row][i];
}
void setequal2(double mat1[n][n], double mat2[n][n], int row)
{
	for (int i = 0; i < n; ++i)
		mat1[row][i] = mat2[i][0];
}
void set_equal(double mat[val][n], double mat1[n][n], int iter)
{
	for (int i = 0; i < n; ++i)
		mat[iter][i] = mat1[0][i];
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

void logic()
{
	for (int i = 0; i < n; ++i)
	{
		u[0][i] = f(x[i]);
		U[0][i] = u[0][i];
	}
	for (int j = 0; j < val; ++j)
	{
		double temp[n][n], temp1[n][n], temp2[n][n];
		setzero(temp);	setzero(temp1);	setzero(temp2);
		setequal(temp, u, 0);
		
		matmult(a2, temp, temp1);
		matadd(temp1, b, temp2);
		setzero(u);
		setequal2(u, temp2, 0);
		set_equal(U, u, j + 1);
	}
}

void euler()
{
	check_gain();//check stability
	boundsandx();//fills boundary, x arrays
	spdiag(); //creates the diagonalized matrix
	logic();
	print(U);
}

int main()
{
	euler();
	getchar();
	return EXIT_SUCCESS;
}

```


### Last Modified:
April 28, 2018
