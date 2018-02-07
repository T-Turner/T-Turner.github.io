---
permalink: /thomasalg
layout: default
---

### Name: thomasalg
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program uses the Thomas Algorithm to solve a tridiagonal matrix and will compute the time needed to solve the matrix. The size of matrix n and time have a linear relationship of time = 2E-8 x n. 

### Input:
The program requires an input for initial conditions and step size.

### Output: 
The program outputs the nxn matrix to time trend to a .csv file and will output the solution to the linear system. 

### Usage:
This program requires that it be called in the following fashion:
Just the Thomas Algorithm:
```c++
  thomas(n);//n is number of rows/columns
```
Time to do the Thomas Algorithm:
```c++
time();
```

### Implementation/Code:
The code is implemented in the following manner:
```c++

#include<iostream>
#include <chrono>
#include<fstream>//file operators

namespace
{
	const int p = 1000;
	double d[p], s[p], l[p], b[p], f[p], x[p];
	double snew[p], bnew[p];
	double h = 0.1, u_a = 5, u_b = 6;
	std::ifstream fin;//file operator
	std::ofstream fout;
}

void setup(int n)
{
	for (int i = 0; i < n; ++i)
	{
		d[i] = -2.0;
		s[i] = 1.0;
		l[i] = 1.0;
	}
	b[0 + 0] = h * h*f[0 + 0] - u_a;
	b[n - 1] = h * h*f[n - 1] - u_b;
	for (int i = 1; i < n-1; ++i)
		b[i] = h * h*f[i];
}
void guassian(int n)
{
	snew[0] = s[0] / b[0];
	bnew[0] = b[0] / d[0];
	for (int i = 1; i < n; ++i)
	{
		double den = d[i] - l[i] * snew[i - 1];
		snew[i] = s[i] / den;
		bnew[i] = (b[i] - l[i] * bnew[i - 1]) / den;
	}
}
void back(int n)
{
	x[n - 1] = bnew[n - 1];
	for (int i = n - 1; i > 0; --i)
		x[i] = bnew[i] - snew[i] * x[i + 1];
}

void print(double mat[], int n)
{
	for (int row = 0; row < n; ++row)
		std::cout << mat[row] << std::endl;
}

void thomas(int n)
{
	setup(n);
	guassian(n);
	back(n);
	//print(x, n);
}

void time()
{
	fout.open("thomas.csv");
	fout << "Number, Time" << std::endl;

	for (int time = 3; time < 1000; ++time)
	{
		auto start = std::chrono::high_resolution_clock::now();
		thomas(time);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		fout << time << "," << elapsed.count() << std::endl;
	}
}

int main()
{
	time();
	std::cout << "FINISHED" << std::endl;
	getchar();
	return EXIT_SUCCESS;
}

```
### Last Modified:
February 7, 2018
