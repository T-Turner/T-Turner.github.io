---
permalink: /implicit_euler
layout: default
---

### Implicit Euler
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a first order ODE.

### Input:
The user must input the left hand side of the equation if the form is y'=a(y)+ect. that is the derivative of y is on the right and the left
side only contains y and t elements. THe user must also input the initial time and final time.

### Output: 
The program will output the solution of the IVP in a table at evenly placed points between the initial and final times

### Usage:

```c++
### Usage:
The constants are defined in the function, and the initial and ending time is defined in the namespace.
```c++
implicit_euler(function name, initial value);
```

### Implementation/Code:
These programs are implemented in the following manner. 

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

	double step = 0.0001;
	double mold, mnew;
	int iter;
	double err;
}

double f1(double x, double y)
{
	return y + dt * (x*lambda) - x;
}
double f2(double x, double y)
{
	return y + dt * (x*gamma - beta * x*x) - x;
}

//NEWTON SOLVER
double deriv(double M, double y, double f(double, double))
{
	return (f(M + step, y) - f(M - step, y)) / (2 * step);
}
void newtlogic(double f(double, double), double y) //newton raphson logic
{
	if (iter != 1)
		mold = mnew;
	if (deriv(mold, y, f) != 0)
		mnew = mold - (f(mold, y) / deriv(mold, y, f));//calculates new guess
	if (mnew != 0 && iter != 1)
		err = abs(mnew - mold) / abs(mnew) * 100; //checks error
}
double newton(double f(double, double), double y)
{
	iter = 0;//resets conditions
	err = 1000;// error is a reset
	mold = 10;
	while (err >= step && iter <= 5000)//continues to do this while () conditions are met
	{
		++iter; //iterates
		newtlogic(f, y);// runs through the logic for the method
	}
	return mnew;
}

//EULER SOLVER
void print1(double yi)
{

	std::cout << "TIME\ty" << std::endl;
	std::cout << ti << "\t" << yi << std::endl;
	for (int i = 1; i <= tf; ++i)
		std::cout << i << "\t" << y[num*i] << std::endl;
}
void implicit_euler(double f(double, double), double yi)
{
	y[0] = yi;
	for (int i = 0; i < n; ++i)
	{
		y[i + 1] = newton(f, y[i]);
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
	implicit_euler(f1, initial_condition);
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
	implicit_euler(f2, initial_condition);
	solver2(g, b, initial_condition);
}

int main()
{
	//comp_exact_euler(1, 5);
	comp_exact_euler2(1, 2, 5);
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
