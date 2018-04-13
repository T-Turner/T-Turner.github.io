---
permalink: /prob6_5
layout: default
---

### Homework 6 Problem #5
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Runs through examples 7.1, 7.2, and 7.3 from "Finite Difference Methods for Ordinary and Partial Differential Equation" 
by Randall J. LeVeque of the University of Washington.

### Input:
No input is required for this program.

### Output: 
Ouputs the Runge Kutta for 3 different ODEs and there respective errors for different time steps

this method is a lot more accurate for all time steps

### Usage:

```c++
example71();
example72();
example73();
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<math.h>
#include<iomanip>

namespace
{
	int n;//iterations
	double yi = 1;//initial conditions
	double y[10000 + 1]{ yi }, k[5][10000 + 1];
	double ti = 0.0, tf = 2;//time interval
	double dt;
}

double f1(double y, double t)
{
	return -sin(t);
}

double f2(double x, double t)
{
	return -10 * (x - cos(t)) - sin(t);
}

double f3(double x, double t)
{
	return -2100 * (x - cos(t)) - sin(t);
}

double fk1(double t, int j, double f(double, double))
{
	return  f(y[j], t);
}

double fk2(double t, int j, double f(double, double))
{
	return f(y[j] + 0.5*k[1][j] * dt, t);
}

double fk3(double t, int j, double f(double, double))
{
	return f((y[j] + 0.5*k[2][j] * dt), t);
}

double fk4(double t, int j, double f(double, double))
{
	return f((y[j] + k[3][j] * dt), t);
}

void rk4(double f(double, double))
{
	for (int j = 0; j < n; ++j)
	{
		k[1][j] = fk1((ti + j * dt), j, f);
		k[2][j] = fk2((ti + j * dt) + dt * 0.5, j, f);
		k[3][j] = fk3((ti + j * dt) + dt * 0.5, j, f);
		k[4][j] = fk4((ti + j * dt) + dt * 1.0, j, f);
		y[j + 1] = y[j] + dt * (k[1][j] / 6 + k[2][j] / 3 + k[3][j] / 3 + k[4][j] / 6);
	}
	std::cout << "RK " << y[n] << std::endl;
}

void example71()
{
	n = 2000;
	dt = (tf - ti) / n;
	rk4(f1);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	std::cout << " Total Error: " << cos(2) - y[n] << std::endl;

}

void example72()
{
	n = 2000;
	dt = (tf - ti) / n;
	rk4(f2);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	std::cout << " Total Error: " << cos(2) - y[n] << std::endl;
}

void example73()
{
	for (int i = 0; i < 5; ++i)
	{
		if (i == 0)
			n = tf / 0.001;
		if (i == 1)
			n = tf / 0.000976;
		if (i == 2)
			n = tf / 0.000950;
		if (i == 3)
			n = tf / 0.0008;
		if (i == 4)
			n = tf / 0.0004;
		dt = (tf - ti) / n;
		rk4(f3);
		std::cout << " Total Error: " << cos(2) - y[n] << std::endl;
	}
}

int main()
{
	std::cout << "EXAMPLE #1" << std::endl;
	example71();
	std::cout << "EXAMPLE #2" << std::endl;
	example72();
	std::cout << "EXAMPLE #3" << std::endl;
	example73();
	getchar();
	return EXIT_SUCCESS;
}
```


### Last Modified:
April 12, 2018
