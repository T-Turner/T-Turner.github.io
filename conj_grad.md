---
permalink: /conj_grad
layout: default
---

### Conjugate Gradient Method
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a sytem of linear equations. It can converge very quickly, but is more unstable than Guass-Seidel or Jacobi iteration.

### Input:
Requires the input of an 'a' and 'b' matrix.
IN the example:
4x+1y=1
x+3y=2
### Output: 
Outputs the solution into the first column of an array 'x'.
The program outputs:
x=1/11
y=7/11

### Usage:

```c++
conj(a, b);
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//Trevor Turner

#include<iostream>
#include<iomanip>

namespace
{
	const int n = 3;
	double a[n][n] = { 10,-5,-4,-5,12,-6,-4,-6,10 };
	double b[n][n] = { 10,0,0,-20,0,0,15,0,0 };
	bool add = true, sub = false;
	double r[n][n], x[n][n], alpha[n][n], trans[n][n], temp[n][n], p[n][n];
	double temp1[n][n], temp2[n][n], rnew[n][n], rold[n][n], beta[n][n];
	double error = 0.000001;
}

void zero(double mat[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			mat[row][col] = 0;
}

void matmult(double mat1[n][n], double mat2[n][n], double result[n][n])
{
	zero(result);
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			for (int last = 0; last < n; ++last)
				result[row][col] += mat1[row][last] * mat2[last][col];
}

void mataddsub(double mat1[n][n], double mat2[n][n], double result[n][n], bool op)
{
	zero(result);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (op == add)
				result[i][j] = mat1[i][j] + mat2[i][j];
			else
				result[i][j] = mat1[i][j] - mat2[i][j];
		}
}

void transpose(double mat1[n][n], double result[n][n])
{
	zero(result);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			result[i][j] = mat1[j][i];
}

void setequal(double blank[n][n], double value[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			blank[i][j] = value[i][j];
}

void print(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			std::cout << mat[i][j] << " ";
		std::cout << std::endl;
	}
}

void f_alpha(double mata[n][n])//works
{
	transpose(p, trans);
	matmult(trans, mata, temp);
	matmult(temp, p, temp2);
	alpha[0][0] = rold[0][0] / temp2[0][0];
	for (int i = 1; i<n; ++i)
		alpha[i][i] = alpha[0][0];
}

void rxnew(double mata[n][n])
{
	//find x
	matmult(alpha, p, temp);
	mataddsub(x, temp, temp1, add);
	setequal(x, temp1);
	//find r
	matmult(alpha, mata, temp);
	matmult(temp, p, temp1);
	mataddsub(r, temp1, temp2, sub);
	setequal(r, temp2);
	transpose(r, trans);
	matmult(trans, r, rnew);
}

void pnew()
{
	beta[0][0] = rnew[0][0] / rold[0][0];//beta
	for (int i = 1; i < n; ++i)
		beta[i][i] = beta[0][0];
	matmult(beta, p, temp1);
	mataddsub(r, temp1, p, add);
	setequal(rold, rnew);
}

void start(double mata[n][n], double matb[n][n])//works
{
	matmult(mata, x, temp);
	mataddsub(matb, temp, r, sub);
	setequal(p, r);
	transpose(r, temp);
	matmult(temp, r, rold);
}

void conj(double mata[n][n], double matb[n][n])
{
	start(mata, matb);
	do
	{
		f_alpha(mata);
		rxnew(mata);
		pnew();
	} while (abs(rnew[0][0]) > error);
	print(x);
}

int main()
{
	conj(a, b);
	getchar();
	return EXIT_SUCCESS;
}
```


### Last Modified:
March 21, 2018
