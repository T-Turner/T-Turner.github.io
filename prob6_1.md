---
permalink: /prob6_1
layout: default
---

### Homework 6 Problem #1
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Runs through examples 7.1, 7.2, and 7.3 from "Finite Difference Methods for Ordinary and Partial Differential Equation" by Randall J. LeVeque of the University of Washington.

### Input:
No input is required for this program.

### Output: 
Ouputs the Euler and Exact solution for 3 different ODEs and there respective errors for different time steps.
The first two examples the output matches the book, but the third example there is some truncation errors which results in different values of the same magnitude.

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
double f2(double x, double y, double t)
{
	return y + dt * (-10 * (y - cos(t)) - sin(t)) - x;
}
double f3(double x, double y, double t)
{
	return y + dt * (-2100 * (y - cos(t)) - sin(t)) - x;
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

void example71()
{
	n = 2000;
	dt = (tf - ti) / n;
	num = n / tf;
	explicit_euler(f1, 1, 2);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	std::cout << " Total Error: " << cos(2) - y[n] << std::endl;

}

void example72()
{
	n = 2000;
	dt = (tf - ti) / n;
	num = n / tf;
	explicit_euler(f2, 1, 2);
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
		num = n / tf;
		explicit_euler(f3, 1, 2);
		std::cout << "Total Error:             " << cos(2) - y[n] << std::endl;
	}
}

int main()
{
	std::cout << "EXAMPLE 7.1" << std::endl;
	example71();
	std::cout << std::endl << "EXAMPLE 7.2" << std::endl;
	example72();
	std::cout << std::endl << "EXAMPLE 7.3" << std::endl;
	example73();
	getchar();
	return EXIT_SUCCESS;
}

```


### Last Modified:
April 12, 2018
