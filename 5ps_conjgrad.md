---
permalink: /5ps_conjgrad
layout: default
---

### 5 Point Stencil Conjugate Gradient
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves the 5 point stencil with the conjugate gradient method.

### Input:
Input boundary conditions.

### Output: 
Solves the 5 point stencil for arbitrary boundary conditions

### Usage:

```c++
conj(A, newB);
```


### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<iomanip>

namespace
{
	const int num = 5; const int m = num - 1;
	const int N = m * m - 2 * m + 2;
	const int n = N - 1;
	double fill[num][num], a[N][N], b[N];
	double A[n][n], B[n], newB[n][n];
	double top = 0, bot = 0,/* rside = 5,*/ lside = 0;
	//conjugate gradient stuff
	bool add = true, sub = false;
	double r[n][n], x[n][n], alpha[n][n], trans[n][n], temp[n][n], p[n][n];
	double temp1[n][n], temp2[n][n], rnew[n][n], rold[n][n], beta[n][n];
	double error = 0.000001;
}
//stencil
void fill_bound(double mat[num][num])//no corners
{
	for (int i = 1; i < m; ++i)
	{
		mat[0][i] = top;
		mat[m][i] = bot;
		mat[i][0] = lside;
		mat[i][m] = sin(5);//rside;
		for (int j = 1; j < m; ++j)
			mat[i][j] = 1;
	}
}
void five_ps(double mata[N][N], double matb[N], double matf[num][num])
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
void clean(double mata[N][N], double matb[N], double newa[n][n], double newb[n])
{
	for (int i = 1; i < N; ++i)
	{
		newb[i - 1] = matb[i];
		for (int j = 1; j < N; ++j)
			newa[i - 1][j - 1] = mata[i][j];
	}
}
void stencil(double mata[N][N], double matb[N], double newa[n][n], double newb[n])
{
	fill_bound(fill);
	five_ps(mata, matb, fill);
	clean(a, b, A, B);
}
//conj grad
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
		std::cout << std::setprecision(6);
		std::cout << mat[i][0] << " ";
		if ((i + 1) % (m - 1) == 0)
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
void adjust(double oldb[n], double newb[n][n])
{
	for (int i = 0; i < n; ++i)
		newb[i][0] = oldb[i];
}
void conj(double mata[n][n], double matb[n][n])
{
	adjust(B, newB);
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
	stencil(a, b, A, B);
	conj(A, newB);

	getchar();
	return EXIT_SUCCESS;
}
```


### Last Modified:
March 21, 2018
