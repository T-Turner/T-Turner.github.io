---
permalink: /prob6_2
layout: default
---

### Homework 6 Problem #2
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Runs through examples 7.1, 7.2, and 7.3 from "Finite Difference Methods for Ordinary and Partial Differential Equation"
by Randall J. LeVeque of the University of Washington.

### Input:
No input is required for this program.

### Output: 
Ouputs the Euler and Exact solution for 3 different ODEs and there respective errors for different time steps.
Unlike the explicit euler the implicit Euler gives much better results, although still large. 

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

//NEWTON SOLVER
double deriv(double M, double y, double f(double, double, double), double t)
{
	return (f(M + step, y, t) - f(M - step, y, t)) / (2 * step);
}
void newtlogic(double f(double, double, double), double y, double t) //newton raphson logic
{
	if (iter != 1)
		mold = mnew;
	if (deriv(mold, y, f, t) != 0)
		mnew = mold - (f(mold, y, t) / deriv(mold, y, f, t));//calculates new guess
	if (mnew != 0 && iter != 1)
		err = abs(mnew - mold) / abs(mnew) * 100; //checks error
}
double newton(double f(double, double, double), double y, double t)
{
	iter = 0;//resets conditions
	err = 1000;// error is a reset
	mold = 10;
	while (err >= step && iter <= 5000)//continues to do this while () conditions are met
	{
		++iter; //iterates
		newtlogic(f, y, t);// runs through the logic for the method
	}
	return mnew;
}

//EULER SOLVER
void implicit_euler(double f(double, double, double), double yi, double time)
{
	y[0] = yi;
	for (int i = 0; i < n + 1; ++i)
		y[i + 1] = newton(f, y[i], i*dt);
	int value = num * time;
	std::cout << "Implicit Euler Solution: " << y[value] << std::endl;
}

void example71()
{
	n = 2000;
	dt = (tf - ti) / n;
	num = n / tf;
	implicit_euler(f1, 1, 2);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	std::cout << " Total Error: " << cos(2) - y[n] << std::endl;

}

void example72()
{
	n = 2000;
	dt = (tf - ti) / n;
	num = n / tf;
	implicit_euler(f2, 1, 2);
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
		implicit_euler(f3, 1, 2);
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
