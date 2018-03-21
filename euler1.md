---
permalink: /euler1 approximation
layout: default
---

### Euler 1st order 
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a 1st order IVP ODE.

### Input:
Input y at time initial, starting and ending time, number of iterations, and equation.

### Output: 
Outputs the y at several time periods

### Usage:

```c++
euler1();
```


### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//y'+ay+b=0

#include<iostream>
#include<math.h>
#include<iomanip>

namespace
{
	const int n = 10000;
	double ti = 0, tf = 10;
	int num = n / tf;
	double yi = -0.4;//intial conditions
	double dt = (tf - ti) / n;
	double y[n];
}

double a(double t)
{
	return 2;
}
double b(double t)
{
	return -3 * cos(t)*exp(-5 * t);
}
void print1()
{
	std::cout << "TIME\ty" << std::endl;
	std::cout << std::setprecision(4);
	std::cout << ti << "\t" << yi << std::endl;
	for (int i = 1; i <= tf; ++i)
		std::cout << i << "\t" << y[num*i] << std::endl;
}

void euler1()
{
	y[0] = yi;
	double t = ti;
	for (int i = 0; i < n; ++i)
	{
		t += dt;
		y[i + 1] = y[i] + dt * (-a(t)*y[i] - b(t));
	}
	print1();
}

int main()
{
	euler1();
	
	getchar();
	return EXIT_SUCCESS;
}
```


### Last Modified:
March 21, 2018
