---
permalink: /epsilon
layout: default
---

## Machine Epsilon

### Name: Machine Epsilon

### Author: Trevor Turner

### Language: C++

### Description/Purpose: 
This program finds the  maximum relative error in a machine, computer, that occurs do to roundoff error. This program only needs to be run once on a given computer as the machine epsilon will not change.

### Input:
There is no required input, but the program analyzes two double type numbers.

### Output: 
The program will return a double that is the realative error in the computer that the program is run on.

### Usage:
This program requires that it be called in the follwowing fashion:
```c++
find_epsilon();	
```
This will automatically print to the screen the machine epsilon.

### Implementation/Code:
This code requires two header files. The two header files allow the function to print to the screen and determine the proper number of significant figures.
```
#include<iostream>//allows input and output
#include<iomanip>//allows for significant figures to be adjusted
```
The code is implemented in the following manner:
```c++
void find_epsilon()//fuction
{
	double epsilon, mach_epsilon = 2;	//initialize values, and sets mach_epsilon to an arbitrary value
	while (1 + mach_epsilon != 1)		//iterates to make sure error is not zero (which is impossible)
	{
		epsilon = mach_epsilon;		//saves previous epsilon
		mach_epsilon /= 2;		//goes to the next step
	}
	std::cout << std::setprecision(6);	//sets the number of significant figures output
	std::cout << "Machine Epsilon is: ";	//output
	std::cout << 2 * epsilon << std::endl;	//output info
}
```

### Last Modified:
January 19, 2018
