---
permalink: /adambash
layout: default
---

### Adams Bashwarth Implicit
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
corrected_adambash(initial value, function);
```


### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//adamds bashworth
#include<iostream>

namespace
{
	const int n = 10000;
	double ti = 0, tf = 10, t;
	double dt = (tf - ti) / n;
	double y[n + 1];
}

double f1(double x, double t)
{
	return -1 * x;
}

double f2(double x, double t)
{
	return 1*x - 5*x*x;
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

void adambash(double x, double f(double, double))
{
	y[0] = x;
	for (int i = 1; i <= n; ++i)
	{
		t = i * dt;
		y[i] = step4(y[i - 1], t, f);
	}
	for (t = ti; t <= tf; t += 1)
	{
		int i = t * n / 10;
		std::cout << "Time " << t << " " << y[i]<< std::endl;
	}
}

void corrected_adambash(double x, double f(double, double))
{
	y[0] = x;
	for (int i = 1; i <= n; ++i)
	{
		t = i * dt;
		y[i] = step4_1(y[i - 1], t, f);
	}
	for (t = ti; t <= tf; t += 1)
	{
		int i = t * n / 40;
		std::cout << "Time " << t << " " << y[i] << std::endl;
	}
}

//u'=lambda*u
void solver1(double l, double alpha)
{
	std::cout << "Exact Solution:" << std::endl;
	for (int time = 0; time <= tf; ++time)
	{
		double result = alpha * exp(time*l);
		std::cout << "Time: " << time << "  " << result << std::endl;
	}
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
int main()
{
	corrected_adambash(5, f1);
	solver1(-1, 5);
	//adambash(25, f2);
	//solver2(1, 5, 25);

	std::cout << "Program has reached the end of life" << std::endl;
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
