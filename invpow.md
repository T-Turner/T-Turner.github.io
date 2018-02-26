---
permalink: /invpow
layout: default
---

## Inverse Power rule
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Finds the smallest and largest eigen vector of a matrix as well as the condition number of the matrix.
The program below uses the 2x2 matrix (1,5,3,3) and outputs eigen vectors of 6, and 2 and a condition number of 3.

### Input:
Input square matrix and number of iterations.

### Output: 
The smallest eigen value and vector and the condition number.

### Usage:
This program requires that it be called in the following fashion:
```c++
findeigen_cond(100);//100 is the number of iterations
```

### Implementation/Code:
The code is implemented in the following manner:
```c++
#include<iostream>
#include<random>

namespace
{
	const int n = 2, m = n - 1;
	double b[n][n], btemp[n][n], bnew[n][n], lambda, save;
	double a[n][n] = { 1,5,3,3 };
	double low[n][n], up[n][n], z[n], y[n];
}

int random()
{
	static std::random_device rd;
	static std::mt19937 mt(rd());
	std::uniform_int_distribution<> dist(1, 9);
	return dist(mt);
}

double norm(double mat[n][n])
{
	double max = 0;
	for (int row = 0; row < n; ++row)
		if (mat[row][0] > max)
			max = mat[row][0];
	return max;
}

void matmult(double mat1[n][n], double mat2[n][n])
{
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			bnew[row][col] = 0;
	for (int row = 0; row < n; ++row)
		for (int col = 0; col < n; ++col)
			for (int last = 0; last < n; ++last)
				bnew[row][col] += mat1[row][last] * mat2[last][col];
}

void setup(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
		mat[i][0] = random();
}

void findeigenvect(int num, double mata[n][n])
{
	setup(b);
	for (int i = 0; i < num; ++i)
	{
		matmult(mata, b);//returns bnew as the new matrix
		lambda = norm(bnew) / norm(b);
		for (int j = 0; j < n; ++j)
			b[j][0] = bnew[j][0];
	}
	std::cout << "Lambda: " << lambda << std::endl;
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

void lu_decom(double mata[n][n], double matb[n][n])
{
	lu_tri(mata);
	find_y(matb);
	find_z();
}

void setequal(double matb[n][n])
{
	for (int i = 0; i < n; ++i)
		matb[i][0] = z[i];
}

void invpow(int num, double mata[n][n])
{
	setup(b);
	for (int i = 0; i < num; ++i)
	{
		lu_decom(a, b);//outputs to z array
		setequal(btemp);//sets z array to bnew
		matmult(a, btemp);//outputs to bnew
		lambda = norm(bnew) / norm(b);
		for (int j = 0; j < n; ++j)
			b[j][0] = bnew[j][0];
	}
	std::cout << "Lambda: " << lambda << std::endl;
}

void findeigen_cond(int num)
{
	findeigenvect(num, a);//outputs a eigen vector in first column b[any integer][0] and eigen value as lambda
	save = lambda;
	invpow(num, a);
	std::cout << "Condition number: " << save / lambda << std::endl;
}


int main()
{
	findeigen_cond(100);
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
February 26, 2018
