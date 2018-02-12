 ---
permalink: /nmnorm
layout: default
---

## Norms
### Name: Non-square matrix norms
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
These programs solve the norms for a given non-square matrix.

### Input:
The program requires input of a mxn matrix.

### Output: 
These program will return a value that is the scalar matrix norm.

### Usage:
P 1 norm:
```c++
p1norm(a);//a is the matrix, naming is arbitrary
```
P  infinity norm:
```c++
pinfnorm(a);a is the matrix, naming is arbitrary
```

### Implementation/Code:
These programs are implemented in the following manner:
double p1norm(double a[m][n])
{
	double max;
	for (int row = 0; row < m; ++row)
	{
		double sum = 0.0;
		for (int col = 0; col < n; ++col)
			sum += a[row][col];
		if (sum > max)
			max = sum;
	}
	return max;
}

double pinfnorm(double a[m][n])
{
	double max;
	for (int col = 0; col < n; ++col)
	{
		double sum = 0.0;
		for (int row = 0; row < m; ++row)
			sum += a[row][col];
		if (sum > max)
			max = sum;
	}
	return max;
}


### Last Modified:
February 12, 2018
