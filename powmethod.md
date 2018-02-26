---
permalink: /powmethod
layout: default
---

## Power Method with Hilbert Matrix
 
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Finds the eigen vector of any square matrix. This program also makes the hilbert matrix.

### Input:
Any square matrix can take the place of the hilbert matrix, it is also required to enter the number of iterations desired.

### Output: 
The program outputs the largest eigen value. for the 3x3 hilbert matrix the eigen value is 1.40832.

### Usage:

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
	const int n = 3, m = n - 1;
	double b[n][n], bnew[n][n], lambda;
	double a[n][n];
	//double a[n][n] = { 1,5,3,3 };
}

void hilbert(double mat[n][n])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			mat[i][j] = pow(i + j + 1, -1);
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

double frobnorm(double mat[n][n])
{
	double fnorm = 0;
	for (int i = 0; i < n; ++i)
		fnorm += pow(mat[i][0], 2);
	return sqrt(fnorm);
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

void print(int iter)//for testing
{
	std::cout << "Iteration: " << iter + 1 << " Lambda: " << lambda << " vector: ";
	for (int j = 0; j < n; ++j)
		std::cout << b[j][0] << " ";
	std::cout << "norm :" << norm(b) << std::endl;
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
	double nnorm = frobnorm(b);
	for (int k = 0; k < n; ++k)
		std::cout << b[k][0] / nnorm << " ";
	std::cout << std::endl << "Lambda: " << lambda << std::endl;
}

int main()
{
	hilbert(a);

	findeigenvect(100);//outputs a eigen vector in first column b[any integer][0] and eigen value as lambda
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
February 26, 2018
