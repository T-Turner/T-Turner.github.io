---
permalink: /ode2
layout: default
---

##Solution of dP/dt=aP+b*P^2

### Name: Ode2

### Author: Trevor Turner

### Language: C++

### Description/Purpose: 
This program finds the solution for the differential equation a range of time steps using the Runge Kutta 4th order method.

### Input:
The program requires the input starting position and velocity of spring mass damper system, also requires values for mass, damping constant, spring constant and the startign and end tiem of the desired interval.

### Output: 
The program will return a range of y and y' values for the given time interval.

### Usage:
This program requires that it be called in the following fashion:
```c++
rk();
```
### Implementation/Code:
The code is implemented in the following manner:
```c++
//y"+2y'+10y=3cos(t)
//y"+ay'+by+c=0
#include<iostream>
#include<math.h>
#include<iomanip>

namespace
{
	const int n = 10000;//iterations
	double yi, ydoti, m, c, k, f;//initial conditions
	int order = 2;
	double y[2][n+1]{ { yi },{ ydoti } }, k1[2][n+1], k2[2][n+1], k3[2][n+1], k4[2][n+1];
	double ti = 0.0, tf = 10.0;//time interval
	double dt = (tf - ti) / n;
}

double coef(char i, double t)
{
	yi = 0;
	ydoti = 0;
	m = 1;
	c = -2;
	k = -10;
	f = 3 * cos(t);
	if (i == 'a')
		return c / m;
	else if (i == 'b')
		return k / m;
	else
		return f / m;
}

double fk1(double t, int i, int j)
{
	if (i == 1)
		return coef('c',t) + coef('a', t) * y[1][j] + coef('b', t) * y[0][j];
	else //i==0
		return y[1][j];
}

double fk2(double t, int i, int j)
{
	if (i == 1)
		return coef('c', t) + coef('a', t) * (y[1][j] + 0.5*k1[1][j] * dt) + coef('b', t) * (y[0][j] + 0.5*k1[0][j] * dt);
	else //i==0
		return y[1][j];
}

double fk3(double t, int i, int j)
{
	if (i == 1)
		return coef('c', t) + coef('a', t) * (y[1][j] + 0.5*k2[1][j] * dt) + coef('b', t) * (y[0][j] + 0.5*k2[0][j] * dt);
	else //i==0
		return y[1][j];
}

double fk4(double t, int i, int j)
{
	if (i == 1)
		return coef('c', t) + coef('a', t) * (y[1][j] + k3[1][j] * dt) + coef('b', t) * (y[0][j] + k3[0][j] * dt);
	else //i==0
		return y[1][j];
}

void rk()
{
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < order; ++i)
		{
			k1[i][j] = fk1((ti + j*dt), i, j);
			k2[i][j] = fk2((ti + j*dt) + dt*0.5, i, j);
			k3[i][j] = fk3((ti + j*dt) + dt*0.5, i, j);
			k4[i][j] = fk4((ti + j*dt) + dt*1.0, i, j);
			y[i][j + 1] = y[i][j] + dt*(k1[i][j] / 6 + k2[i][j] / 3 + k3[i][j] / 3 + k4[i][j] / 6);
		}
}

void print()
{
	std::cout << "Runge Cutta Method for 2nd Order ODE:" << std::endl << std::endl;
	std::cout << "TIME\tY   \t\tYdot   " << std::endl;
	std::cout << std::setprecision(5);
	for (int j = 0; j <= n; j += n / 10)
		std::cout << ti + j*dt << "\t" << y[0][j] << "\t" << y[1][j] << std::endl;
}

int main()
{
	rk();
	print();
	getchar();
	return EXIT_SUCCESS;
}
```

### Last Modified:
January 22, 2018
