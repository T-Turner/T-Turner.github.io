---
permalink: /guass_seidel
layout: default
---

### Guass-Seidel
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a diagonally dominant, or symmetric and positive definite system of equations.

### Input:
requires the user to input the 'A' and 'B' matricies. 

### Output: 
Outputs the solution of the linear system into a 1-D array called ans.

### Usage:

```c++
  guasssied()
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//Trevor Turner
//A01975883
//guass siedel's convergence is only guaranteed if the matrix is either diagonally dominant, or symmetric and positive definite.
#include<iostream>
#include<iomanip>

namespace
{
	const int n = 3;
	double a[n][n] = { 3,-0.1,-0.2, 0.1, 7.0, -0.3,0.3,-0.2,10 };
	double b[n] = { 7.85, -19.3, 71.4 };
	double x[n], ans[n];
	double err[n], sum=100;
	double error = 0.0001;
}

void populate()
{
	for (int i = 0; i < n; ++i)
		err[i] = 100;
}

void check(int i)
{
	x[i] = ans[i];
	ans[i] = b[i];
	for (int j = 0; j < i; ++j)
		ans[i] += -a[i][j] * ans[j];
	for (int j = i + 1; j < n; ++j)
		ans[i] += -a[i][j] * ans[j];
	ans[i] = ans[i] / a[i][i];
	err[i] = abs((ans[i] - x[i]) / ans[i]);
}

void print()
{
	for (int i = 0; i < n; ++i)
		std::cout << ans[i] << " ";
	std::cout << std::endl;
}

void guasssied()
{
	populate();
  while ((sum/n) > error)
	{
		sum = 0;
		for (int i = 0; i < n; ++i) 
		{
			sum += err[i];
			if (err[i] > error)
				check(i);
		}
	}
  print();
}

int main()
{
	guasssied();
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
March 21, 2018
