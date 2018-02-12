 ---
permalink: /norms
layout: default
---

## Norms
### Name: P norms and Frobinius norm for nxn matrix
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
These programs solve the norms for a given matrix.

### Input:
THe p norm requires one to enter which p norm and the 'a' matrix, the frobinius norm requires just the 'a' matrix.

### Output: 
These program will return the norm under the variable name of 'max'.

### Usage:
The p norms are called:
```c++
pnorm(p, a);//"p" is which norm and a is the matrix
```
THe frobinius norm is called:
```c++
frobnorm(a);//'a' is the matrix
```

### Implementation/Code:
These programs are implemented in the following manner:
P norms:
```c++
namespace
{
	const int n = 10;
	double norm[n];
	double pi = 3.1415926535;
}

double pnorm(int p, double a[n][n])
{
	double sum, max = 0.0;
	if (p == 1)
		p1norm(sum, max, a);
	else if(p ==2)
		p2norm(sum, max, a);
	else
		pinfnorm(sum, max, a);
	return max;
}

void p1norm( double sum, double max, double a[n][n])
{
	for (int row = 0; row < n; ++row)
	{
		sum = 0.0;
		for (int col = 0; col < n; ++col)
			sum += a[row][col];
		if (sum > max)
			max = sum;
	}
}

void p2norm(double sum, double max, double A[n][n])
{
	double a = A[0][0];
	double b = A[0][1];
	double c = A[1][0];
	max = 2 * sqrt(b*c)*cos(pi / (n + 1)) + a;
}

void pinfnorm(double sum, double max, double a[n][n])
{
	for (int col = 0; col < n; ++col)
	{
		sum = 0.0;
		for (int row = 0; row < n; ++row)
			sum += a[row][col];
		if (sum > max)
			max = sum;
	}
}
```
frobinius norm:
```c++

double frobnorm(double a[n][n])
{
	double norm;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			norm += pow(a[i][j], 2);
	return sqrt(norm);
}
```

### Last Modified:
February 12, 2018
