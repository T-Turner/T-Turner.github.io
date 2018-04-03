---
permalink: /explicit_euler
layout: default
---

### Explicit Euler
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a first order ODE.

### Input:
The user must input the left hand side of the equation if the form is y'=a*y+ect. that is the derivative of y is on the right and the left
side only contains y and t elements. THe user must also input the initial time and final time.

### Output: 
The program will output the solution of the IVP in a table at evenky placed points between the initial and final times

### Usage:
The constants are defined in the function, and the initial and ending time is defined in the namespace.
```c++
explicit_euler(function name, initial value);
```

### Implementation/Code:
These programs are implemented in the following manner. In this form, it compares the exact value to the euler solution.

```c++
#include<iostream>
#include<iomanip>

namespace
{
	const int n = 10000;
	double ti = 0, tf = 10;
	int num = n / tf;
	double dt = (tf - ti) / n;
	double y[n], lambda, beta, gamma;
}

double f1(double y)
{
	return y + dt * (y*lambda);
}
double f2(double y)
{
	return y + dt * (y*gamma - beta * y*y);
}

//EULER SOLVER
void print1(double yi)
{

	std::cout << "TIME\ty" << std::endl;
	std::cout << ti << "\t" << yi << std::endl;
	for (int i = 1; i <= tf; ++i)
		std::cout << i << "\t" << y[num*i] << std::endl;
}
void explicit_euler(double f(double), double yi)
{
	y[0] = yi;
	for (int i = 0; i < n; ++i)
	{
		y[i + 1] = f(y[i]);
	}
	std::cout << "Euler Solution:" << std::endl;
	print1(yi);
}

//u'=lambda*u
void solver1(double l, double alpha)
{
	lambda = l;
	std::cout << "Exact Solution:" << std::endl;
	for (int time = 0; time <= tf; ++time)
	{
		double result = alpha * exp(time*lambda);
		std::cout << "Time: " << time << "  " << result << std::endl;
	}
}
void comp_exact_euler(double l, double initial_condition)
{
	lambda = l;
	explicit_euler(f1, initial_condition);
	solver1(l, initial_condition);
}

//p'=p(gamma-beta)
void solver2(double g, double b, double initial_condition)
{
	std::cout << "Exact Solution:" << std::endl;
	for (int time = 0; time <= tf; ++time)
	{
		double expo = exp(time*g);
		double num = -initial_condition * g*expo;
		double den = initial_condition * b - g - initial_condition * b*expo;
		double result = num / den;
		std::cout << "Time: " << time << "  " << result << std::endl;
	}
}
void comp_exact_euler2(double g, double b, double initial_condition)
{
	beta = b;
	gamma = g;
	explicit_euler(f2, initial_condition);
	solver2(g, b, initial_condition);
}

int main()
{
	comp_exact_euler(-1, 5);
	//comp_exact_euler2(1, 2, 5);
	getchar();
	return EXIT_SUCCESS;
}

```
### Example
For u'=lambda(u) where u(0)=5

lambda=-100, the value of u quickly drops to zero, at 1 second the value is <<<1

lambda=-1, the value of u drops slowly towards zero, at 1 second the value is 1.83

lambda=1, the value of u slowly increases to infinity, at 1 second the value is 13.6

for p'=gamma(p)-beta(p^2); gamma is 0.1, and beta is 0.0001

P(0)=25, the value of u at time=1 is 27.6 and at time=10 65.2

P(0)=40,000, the value of u at time=1 is approx. 8480 and at time=10 approx 1560

### Last Modified:
April 2, 2018
