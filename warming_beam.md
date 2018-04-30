---
permalink: /warming_beam
layout: default
---

### Hypepbolic Equation Beam and Warming
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solves a U_t+2*U_x=0 equation.

### Input:
Requires the initial conditions, the parameters for x and t.

### Output: 
Outputs the solution at a particular point in time and position.

### Usage:

```c++
  beamw();
```


### Implementation/Code:
These programs are implemented in the following manner. 

```c++
#include<iostream>

//ut=-2ux

namespace
{
	const double pi = 3.1415926535;
	const int n = 5;
	const int time = 100;
	double u[time][n + 1];

	double a = -2;
	double x0 = 0, xf = 1;
	double t0 = 0, tf = 2;

	double dx = (xf - x0) / n;
	double dt = (tf - t0) / time;
	double c = a*dt / dx;
}

void stabil_check()
{
	if (abs(c) > 1)
		std::cout << "UNSTABLE TIME STEP; THE RESULTS MAY BE INACCURATE" << std::endl;
}
void bounds()
{
	for (int i = 0; i < time; ++i)
	{
		//u[0][i] = 0;
		u[i][0] = 5;
		u[i][n] = 5;
	}
}
void print()
{
	for (int i = 0; i < time; i += time / 5)
	{
		for (int j = 0; j < n + 1; ++j)
			std::cout << u[i][j] << " ";
		std::cout << std::endl;
	}
}

void beamw()
{
	for (int i = 1; i < time; ++i)
		for (int j = 1; j < n; ++j)
			u[i][j] = u[i - 1][j] - (c / 2)*(3*u[i - 1][j] - 4*u[i - 1][j - 1]+u[i-1][j-2])
			+ (.5*c*c)*(u[i - 1][j] - 2 * u[i - 1][j - 1] + u[i - 1][j - 2]);
	print();
}

int main()
{
	stabil_check();
	bounds();
	beamw();
	getchar();
	return EXIT_SUCCESS;
}

```

### Last Modified:
April 30, 2018
