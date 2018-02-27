---
permalink: /eigen4fd
layout: default
---

## Eigen value for 2nd order Finite Difference
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Finds the eigne value of the 2nd order finite difference matrix. 3=<Lambda<=4
2x2 matrix lambda=3
3x3 matrix lambda=2+sqrt(2)
4x4 matrix lambda=(5+sqrt(5))/2
5x5 matrix lambda=2+sqrt(3)

### Input:
Basically just input number of iterations and the size of the matrix, matrix A. 

### Output: 
Outputs matrix a if desired and the eigen value, lambda. 

### Usage:
This program requires that it be called in the following fashion:
```c++
findeigenvect(100);
```
### Implementation/Code:
The code is implemented in the following manner:
```c++
#include<iostream>
#include<random>

namespace
{
	const int n = 100, m = n - 1;
	double b[n][n], bnew[n][n], lambda;
	double a[n][n];
}

int random()
{
	static std::random_device rd;
	static std::mt19937 mt(rd());
	std::uniform_int_distribution<> dist(0, 10);
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

void findeigenvect(int num)
{
	setup(b);
	for (int i = 0; i < num; ++i)
	{
		matmult(a, b);//returns bnew as the new matrix
		lambda = norm(bnew) / norm(b);
		for (int j = 0; j < n; ++j)
			b[j][0] = bnew[j][0];
	}
	std::cout << "Lambda: " << lambda << std::endl;
}

void set(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
	{
		mat[i][i] = 2;
		if (i < m)
			mat[i][i + 1] = -1.0;
		if (i > 0)
			mat[i][i - 1] = -1.0;
	}
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
	set(a);
	//print(a);
	findeigenvect(100);//outputs a eigen vector in first column b[any integer][0] and eigen value as lambda
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
february 26, 2018
