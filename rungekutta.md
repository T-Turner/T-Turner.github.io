---
permalink: /rungekutta
layout: default
---

### Runge Kutta 2nd and 4th order
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves an ODE.

### Input:
The left hand side of the function (y'=a*y +etc). Also in the namespace define starting and ending times as well as the initial value.

### Output: 
The program will output the solution of the IVP in a table at evenly placed points between the initial and final times

### Usage:
2nd order Runge Kutta
```c++
rk2(function);
```
4th order Runge Kutta
```c++
rk4(function);
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>
#include<math.h>
#include<iomanip>

namespace
{
	const int n = 10000;//iterations
	double yi = 5;//initial conditions
	double y[n + 1]{ yi }, k[5][n + 1];
	double ti = 0.0, tf = 10.0;//time interval
	double dt = (tf - ti) / n;
}

double f1(double y, double t)
{
	return 1 * y;
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

void rk2(double f(double, double))
{
	for (int j = 0; j < n; ++j)
	{
		k[1][j] = fk1((ti + j * dt), j, f);
		k[2][j] = fk2((ti + j * dt) + dt * 0.5, j, f);
		y[j + 1] = y[j] + .5*dt * (k[1][j] + k[2][j]);
	}
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
}

void print()
{
	std::cout << "Runge Cutta Method for 2nd Order ODE:" << std::endl << std::endl;
	std::cout << "TIME\tY" << std::endl;
	std::cout << std::setprecision(5);
	for (int j = 0; j <= n; j += n / 10)
		std::cout << ti + j * dt << "\t" << y[j] << std::endl;
}

int main()
{
	//rk2(f1);
	rk4(f1);
	print();
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
