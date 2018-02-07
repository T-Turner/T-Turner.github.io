---
permalink: /initialcond
layout: default
---

## Initial Conditions

### Name: initialcond
### Author: Trevor Turner
### Language: C++

### Description/Purpose: 
This program uses a basic approach to acquire initial conditions.

### Input:
The program requires the input to be inserted via the cmd bar.

### Output: 
This program outputs instructions to the cmd bar.

### Usage:
These programs follow the same logic pattern. The relative error function is called via:
```c++
  process_input()
```

### Implementation/Code:
These programs are implemented in the following manner. 
```c++
#include<iostream>
#include<iomanip>

namespace
{
	const int n = 1000;
	double u[n];
	double a, b, u_a, u_b, step, time;
	char f;
}

void input()
{
	std::cout << "Elliptic ODE Solver boundary conditions u(a)=u_a and u(b)=u_b " << std::endl;
	std::cout << std::endl << "Enter condition 'a' time step ";
	std::cin >> a;
	std::cout << std::endl << "Enter 'U_a' ";
	std::cin >> u_a;
	std::cout << std::endl << "Enter condition 'b' time step ";
	std::cin >> b;
	std::cout << std::endl << "Enter 'U_a' ";
	std::cin >> u_b;
	std::cout << "Enter forcing function: cos for cosine, sin for sine, tan for tangent";
	std::cout << "exp for exponential ";
	std::cin >> f;
}

double force(double x)
{
	if (f == 'c')
		return cos(x);
	else if (f == 's')
		return sin(x);
	else if (f == 't')
		return tan(x);
	else if (f == 'e')
		return exp(x);
}

void process_input()
{
	input();
	u[0] = u_a;
	u[n - 1] = u_b;
	step = abs(b - a) / n;
	force(f);
}
```

### Last Modified:
January 19, 2018
