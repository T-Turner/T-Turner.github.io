---
permalink: /ode1
layout: default
---

##Solution of dP/dt=aP+b*P^2

### Name: Ode1

### Author: Trevor Turner

### Language: C++

### Description/Purpose: 
This program finds the solution for the differential equation at any given time step.

### Input:
The program requires the input of the the instance of time you want the solution, variable 'a', variable 'b', and the initial conditon of the ordinary differential equation.

### Output: 
The program will return a double that is the solution at that time step of the ODE.

### Usage:
This program requires that it be called in the follwowing fashion:
```c++
std::cout<<P(double t, double a, double b, double p_0)<<std::endl;
```
This is written as such because it is a double type function.

### Implementation/Code:
The code is implemented in the following manner:
```c++
double P(double t, double a, double b, double p_0)
{
	double constant = log(p_0 / (b*p_0 - a));
	double exponent = -a*(t + constant/a);
	return a / (b - exp(exponent));
}
```

### Last Modified:
January 19, 2018
