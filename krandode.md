---
permalink: /krandode
layout: default
---

## Elliptic ODE Solver
### Name: Elliptic ODE Solver
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program solves the elliptic ode of form d/dxk(x)du/dx=f(x) using a random number generator for k. 

### Input:
Input the intial conditions for the ODE and forcing function.

### Output: 
This program outputs the solution to an elliptic ODE.

### Usage:
This program is called via:
```c++
thomas();
```

### Implementation/Code:
```c++

#include<iostream>
#include<random>

namespace
{
	const int n = 10;
	double d[n], s[n], l[n], b[n], f[n], x[n];
	double snew[n], bnew[n];
	double u_a = 5, u_b = 6;//input initial conditons here
	double h = abs(u_a - u_b) / n;
}

double force(double x)
{
	return x;//input function here!!!!!
}

int k()
{
	static std::random_device rd;
	static std::mt19937 mt(rd());
	std::uniform_int_distribution<> dist(10, 50);
	return dist(mt);
}
double find_f()
{
	for (int i = 0; i < n; ++i)
	{
		f[i] = force(x[i]) / k();
	}
}
void setup()
{
	for (int i = 0; i < n; ++i)
	{
		d[i] = -2.0;
		s[i] = 1.0;
		l[i] = 1.0;
	}
	b[0 + 0] = h * h*f[0 + 0] - u_a;
	b[n - 1] = h * h*f[n - 1] - u_b;
	for (int i = 1; i < n - 1; ++i)
		b[i] = h * h*f[i];
}
void guassian()
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
void back()
{
	x[n - 1] = bnew[n - 1];
	for (int i = n - 1; i > 0; --i)
		x[i] = bnew[i] - snew[i] * x[i + 1];
}
void print(double mat[])
{
	for (int row = 0; row < n; ++row)
		std::cout << mat[row] << std::endl;
}
void thomas()
{
	find_f();
	setup();
	guassian();
	back();
	print(x);
}

int main()
{
	thomas();
	getchar();
	return EXIT_SUCCESS;
}

```


### Last Modified:
February 7, 2018
