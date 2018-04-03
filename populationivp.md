---
permalink: /populationivp
layout: default
---

### First order IVP ODE and Logistic Population model
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
Solver 1 solves the basic first order ODE of the form u'= constant(u).
Solver 2 solves the logistic model for population growth of the form p'=constant1(p)-constant2(p^2)

### Input:
Solver 1 requires the user input the constant, initial condition, and what time one wishes to solve for.
Solver 2 requires the user input the two constants, the initial conditions and and the time one wishes to solve for.

### Output: 
Both solvers will output the solution to their respective ODES.

### Usage:
Solver 1:
```c++
solver1(constant, initial value, time);
```
Solver 2:
```c++
solver2(constant1, constant2, initial value, time);
```

### Implementation/Code:
These programs are implemented in the following manner. 

```c++
//u'=lambda*u
void solver1(double l, double alpha, double time)
{
	std::cout << "solution to IVP: " << std::endl << "Time: " << time << "  " << alpha * exp(time*l) << std::endl;
}

//p'=p(gamma-beta*p)
void solver2(double g, double b, double initial_condition, double time)
{
	double expo = exp(time*g);
	double num = -initial_condition * g*expo;
	double den = initial_condition * b - g - initial_condition * b*expo;
	std::cout << "Solution to Population Logistic Model: " << std::endl << "Time: " << time << "  " << num / den << std::endl;
}
```
### Example

u’=lambda(u)
u(0)=5; lambda=1

(time, value)

(0.0, 5.0), (1.0, 13.6), (2.0, 36.9), (3.0, 100.4), (6.0, 2017.1), (10.0, 110132)

p’=alpha(p)-beta(p^2)
alpha=1; gamma=2; p=5;

(time, value)

(0.0, 5.0), (1.0, 0.75), (2.0, 0.57), (3.0, 0.52), (6.0, 0.51), (10.0, 0.500)

### Last Modified:
April 2, 2018
