---
permalink: /jacobi
layout: default
---

## Jacobi Iteration
### Name: initialcond
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program solves a system of tridiagonal linear equations, and the time will follow a logarithmic path to inifinite as the number of equations increases.

### Input:
The program requires the matrix A, and B be defined in the namespace.
 
### Output: 
This program outputs the solution in a nxn array xnew. The first column is the solution.

### Usage:
```c++
  jaco_iter(double a[n][n], double b[n][n]);
````
### Implementation/Code:
These programs are implemented in the following manner. 
```c++
//Trevor Turner
// Jacobi iteration

#include<iostream>
#include<iomanip>

namespace
{
	const int n = 2;
	double u[n][n], l[n][n], d[n][n], xold[n][n], xnew[n][n];
	bool add = true, sub = false;
	double error = 0.5, err = 1000;
	double A[n][n];//enter diagonally dominant 
	double B[n][n];// enter conditions in the first column
}

void back_sub(double mata[n][n], double matb[n][n], double mat[n][n])
{
	for (int row = n - 1; row >= 0; --row)
	{
		mat[row][n - 1] = matb[row][n - 1];
		for (int col = row + 1; col < n; ++col)
			mat[row][n - 1] += -mata[row][col] * mat[col][n - 1];
		mat[row][n - 1] = mat[row][n - 1] / mata[row][row];
	}
}

void for_sub(double mata[n][n], double matb[n][n], double mat[n][n])
{
	for (int row = 0; row < n; ++row)
	{
		mat[row][n - 1] = matb[row][n - 1];
		for (int col = 0; col < row; ++col)
			mat[row][n - 1] += -mata[row][col] * mat[col][n - 1];
		mat[row][n - 1] = mat[row][n - 1] / mata[row][row];
	}
}

void setequal(double mat1[n][n], double mat2[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			mat1[row][col] = mat2[row][col];
}
void setnum(double mat[n][n], double value)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = value;
}
void inv(double mat[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
		{
			if (mat[row][col] == 0)
				mat[row][col] = 0;
			else
				mat[row][col] = pow(mat[row][col], -1);
		}
}
void lud(double mat[n][n])
{
	for (int row = 0; row < n; ++row)
	{
		for (int col = row + 1; col < n; ++col)
			u[row][col] = mat[row][col];
		for (int col = 0; col < row; ++col)
			l[row][col] = mat[row][col];
		d[row][row] = mat[row][row];
	}
}
void matprint(double mat[n][n])
{
	for (int row = 0; row < n; ++row)
	{
		for (int col = 0; col < n; ++col)
			std::cout << mat[row][col] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
void mataddsub(double mat1[n][n], double mat2[n][n], double result[n][n], bool op)
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
		{
			if (op == add)
				result[row][col] = mat1[row][col] + mat2[row][col];
			else //op==subtract
				result[row][col] = mat1[row][col] - mat2[row][col];
		}
}
void matmult(double mat1[n][n], double mat2[n][n], double result[n][n], int row1, int col1, int row2, int col2)
{
	for (int row = 0; row < row1; ++row)
		for (int col = 0; col < col2; ++col)
			for (int last = 0; last < col1; ++last)
				result[row][col] += mat1[row][last] * mat2[last][col];
}
void jaco_logic(double a[n][n], double b[n][n], double xtemp[n][n])
{
	double xtemp1[n][n], xtemp2[n][n];
	setnum(xtemp1, 0.0);
	setnum(xtemp2, 0.0);
	setequal(xold, xnew);
	matmult(xtemp, xold, xtemp1, n, n, n, 1);
	mataddsub(b, xtemp1, xtemp2, sub);
	inv(d);
	matmult(d, xtemp2, xnew, n, n, n, 1);
}

void jaco_iter(double a[n][n], double b[n][n])
{
	double xtemp[n][n];
	lud(a);
	mataddsub(l, u, xtemp, add);
	while (error<err)
	{
		jaco_logic(a, b, xtemp);
		err = abs(xnew[0][0] - xold[0][0]);
	} 
	matprint(xnew);
	matprint(xold);
}

int main()
{
	jaco_iter(A, B);
	getchar();
	return EXIT_SUCCESS;
}
```
### Last Modified:
February 7, 2018
