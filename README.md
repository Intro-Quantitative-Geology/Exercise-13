# Lab-exercise-6
This exercise is part 1 of the exercises on thermochronology.
In this exercise you will run a transient 1D thermal model, plot predicted geotherms and predict thermochronometer ages for several thermochronometers.

## Overview
This two-part set of exercises is designed to give you a better understanding of thermochronology and the thermal field in the Earth's crust.
Thermochronology combines many aspects of what has been studied in this course, including the advection and diffusion equations, tectonic and surface (erosional) processes and some basic geostatistics.
The first part of the laboratory exercises on thermochronology is intended to provide an understanding of heat advection and diffusion in the Earth's crust, the time-dependence of thermal processes, and how a thermal history can be used to predict thermochronometer ages.

## Time-dependent temperature in the Earth
In this exercise we will use an analytical solution to the 1-D time-dependent thermal advection-diffusion equation to simulate erosion of rock at the Earth's surface, the upward transport of the underlying rock toward the surface and the changes in a 1-D geotherm with time.
The basic equation for temperature *T* as a function of depth *z* and time *t* was originally published by Carslaw and Jaeger (1959)

![Equation 1](Images/Equation1.png)<br/> *Equation 1. 1D time-dependent heat advection-diffusion equation.*

where *G* is the initial geothermal gradient (increase in temperature with depth), *v*<sub>*z*</sub> is the vertical advection velocity (positive upward), *κ* is the thermal diffusivity and `erfc()` is the complementary error function, defined as

![Equation 2](Images/Equation2.png)<br/> *Equation 2. The complementary error function.*

With this equation, the temperature initially increases linearly with depth (*T*(*z*,*t*=0) = *Gz*) and the geotherm will evolve with time as a function of the advection velocity and thermal diffusivity.

## Problem 1 - The time dependence of rock advection
