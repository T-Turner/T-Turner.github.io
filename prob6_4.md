---
permalink: /prob6_4
layout: default
---

### Homework 6 Problem #4
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Runs through examples 7.1, 7.2, and 7.3 from "Finite Difference Methods for Ordinary and Partial Differential Equation" by 
Randall J. LeVeque of the University of Washington.

### Input:
No input is required for this program.

### Output: 
Ouputs the Adam's Bashforth Corrected 4th order for 3 different ODEs and there respective errors for different time steps

unfortunately the multistep method is completely inaccurate...

### Usage:

```c++
example71();
example72();
example73();
```


### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//adams bashforth
#include<iostream>

namespace
{
	double ti = 0, tf = 2, t;
	int n;
	double dt;
	double y[10000 + 1];
}

double f1(double x, double t)
{
	return -sin(t);
}
double f2(double x, double t)
{
	return (-10 * (-cos(t)) - sin(t)) - x;
}
double f3(double x, double t)
{
	return (-2100 * (-cos(t)) - sin(t)) - x;
}

double step1(double x, double t, double f(double, double))
{
	return x + dt * f(x, t);
}
double step2(double x, double t, double f(double, double))
{
	t += dt;
	return step1(x, t, f) + (1.5*f(step1(x, t, f), t) - .5*f(x, t))*dt;
}
double step3(double x, double t, double f(double, double))
{
	t += dt;
	return step2(x, t, f) + ((23 / 12)*step2(x, t, f) - (4 / 3)*step1(x, t, f) + (5 / 12)*f(x, t))*dt;
}
double step4(double x, double t, double f(double, double))
{
	t += dt;
	return step3(x, t, f) + ((55 / 24)*step3(x, t, f) - (59 / 24)*step2(x, t, f) + (37 / 24)*step1(x, t, f) - (3 / 8)*f(x, t))*dt;
}

double step1_1(double x, double t, double f(double, double))
{
	return x + (f(step1(x, t, f), t) + f(x, t))*dt / 2;
}
double step2_1(double x, double t, double f(double, double))
{
	return step1_1(x, t, f) + dt * (5 * f(step2(x, t, f), t) + 8 * f(step1_1(x, t, f), t) - f(x, t)) / 12;
}
double step3_1(double x, double t, double f(double, double))
{
	return step2_1(x, t, f) + dt * (9 * f(step3(x, t, f), t) + 19 * f(step2_1(x, t, f), t) - 5 * f(step1_1(x, t, f), t) + f(x, t)) / 24;
}
double step4_1(double x, double t, double f(double, double))
{
	return step3_1(x, t, f) + dt * (251 * f(step4(x, t, f), t)
		+ 646 * f(step3_1(x, t, f), t) - 264 * f(step2_1(x, t, f), t)
		+ 106 * f(step1_1(x, t, f), t) - 19 * f(x, t)) / 720;
}


void corrected_adambash(double x, double f(double, double))
{
	y[0] = x;
	for (int i = 1; i <= n; ++i)
	{
		t = i * dt;
		y[i] = step4_1(y[i - 1], t, f);
	}
	int value = tf * n /(4*tf);
	std::cout << "Time " << t << " solution: " << y[value] << std::endl;
}

void example71()
{
	n = 2000;
	dt = (tf - ti) / n;
	corrected_adambash(1, f1);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	int value = tf * n / (4*tf);
	std::cout << " Total Error: " << cos(2) - y[value] << std::endl;

}

void example72()
{
	n = 2000;
	dt = (tf - ti) / n;
	corrected_adambash(1, f2);
	std::cout << "Exact          Solution: " << cos(2) << std::endl;
	int value = tf * n / (4 * tf);
	std::cout << " Total Error: " << cos(2) - y[value] << std::endl;
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
		corrected_adambash(1, f3);
		int value = tf * n / (4 * tf);
		std::cout << " Total Error: " << cos(2) - y[value] << std::endl;
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

	std::cout << "Program has reached the end of life" << std::endl;
	getchar();
	return EXIT_SUCCESS;
}
```


### Last Modified:
April 12, 2018
