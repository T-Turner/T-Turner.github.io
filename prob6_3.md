---
permalink: /prob6_3
layout: default
---

### Homework 6 Problem #3
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Runs through examples 7.1, 7.2, and 7.3 from "Finite Difference Methods for Ordinary and Partial Differential Equation"
by Randall J. LeVeque of the University of Washington.

### Input:
No input is required for this program.

### Output: 
Ouputs the explict Euler for the ODE dy/dt=-sin(t) and the respective errors for different time steps.
STEP SIZE: 5000 the error was -0.000181878

STEP SIZE: 2000 the error is -0.000454767

STEP SIZE: 1000 the error is -0.000909769

STEP SIZE: 800 the error is -0.00113736

STEP SIZE: 600 the error is -0.00151681

STEP SIZE: 400 the error is -0.00227619

STEP SIZE: 200 the error is -0.00455829

### Usage:

```c++
stepper();
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<iomanip>

namespace
{
	double ti = 0, tf = 2;
	int n, num;
	double dt;
	double y[10000 + 1];

	double step = 0.0001;
	double mold, mnew;
	double err;
	int iter;
}

double f1(double x, double y, double t)
{
	return y - dt * sin(t) - x;
}

//EULER SOLVER
void explicit_euler(double f(double, double, double), double yi, double time)
{
	y[0] = yi;
	for (int i = 0; i < n + 1; ++i)
		y[i + 1] = f(0, y[i], i*dt);
	int value = num * time;
	std::cout << "Explicit Euler Solution: " << y[value] << std::endl;
}

void logic()
{
	dt = (tf - ti) / n;
	num = n / tf;
	explicit_euler(f1, 1, 2);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	std::cout << " Total Error: " << cos(2) - y[n] << std::endl;
}

void stepper()
{
	n = 5 * 1000;
	std::cout << "STEP SIZE: " << n << std::endl;
	logic();
	n = 2 * 1000;
	std::cout << "STEP SIZE: " << n << std::endl;
	logic();
	
	for (double i = 1; i > 0.1; i -= 0.2)
	{
		n = i * 1000;
		std::cout << "STEP SIZE: " << n << std::endl;
		logic();
	}
}

int main()
{
	stepper();
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
April 12, 2018
