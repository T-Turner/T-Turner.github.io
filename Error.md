---
permalink: /error
layout: default
---

## Absolute and Relative Error

### Name: Error

### Author: Trevor Turner

### Language: C++

### Description/Purpose: 
These programs find the relative and absolute error between a true or theoretical value and an experimental value.

### Input:
Both of these functions require 2 double type numbers, namely the theoretical value and the estimated value.

### Output: 
These program will return a double that is the relative error, or the absolute error.

### Usage:
These programs follow the same logic pattern. The relative error function is called via:
```c++
double rel_error(double actual, double theo)
```
The absolute error is called via:
```c++
double abs_error(double actual, double theo)
```
These functions return a double type such that the error could be analyzed to achieve a certain accuracy. These functions require 2 double type the first is the actual value and the second is the theoretical value.

### Implementation/Code:
These programs are implemented in the following manner. 
The relative error function:
```c++
double rel_error(double actual, double theo)
{
	return abs_error(actual, theo) / abs(actual);
}
```
The absolute error:
```c++
double abs_error(double actual, double theo)
{
	return abs(theo - actual);
}
```

### Last Modified:
January 19, 2018
