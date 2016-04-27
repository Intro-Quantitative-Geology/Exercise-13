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
We will begin by plotting geotherms to get a sense of how our temperature equation works.

1. If you [download a copy of the Python script `age_predict_1D.py`](age_predict_1D.py) you should be able to run it without making any changes to produce a plot like that shown below.

    ![Time-dependent heat transfer with advection](Images/1D_transient_plot1.png)<br/>
    *Figure 1. 1D transient thermal solution including advection.*
To start, please add axis labels and a title to this plot. What is the advection velocity for this thermal solution? Add a text label listing the advection velocity on your the plot using the text() function, then save a copy of the plot for your write-up.
2. How does the thermal solution change when you alter the advection velocity? Increase the advection velocity to 1.0 mm/a and save the resulting plot for your write-up. Do the same thing for an advection velocity of 0.1 mm/a. What is the effect of changing the advection velocity on temperatures in the shallow crust (<10 km depth)?
3. Reset the advection velocity to the starting value. Now increase the total simulation time to 100 Ma. Once again, save a copy of this plot for your write-up. What happens to temperatures in the model as the simulation time increases? You might want to run some additional calculations (you don't need to save the plots) with even longer simulation times (1000 Ma, 5000 Ma, etc.). What do you observe for the temperatures in the model? Are there any potential problems with these temperatures? Does the model approach a steady-state thermally?
4. Reset the simulation time to 50 Ma. Now we'll explore the effect of the initial thermal gradient. Increase the initial thermal gradient to 20°C/km, run the model and save the plot. Do the same for a thermal gradient of 5°C/km. How does the thermal gradient affect the temperatures in the model? Is there any clear relationship between the initial thermal gradient and the maximum temperature in the model at t = 0 Ma?
5. Lastly, reset the initial thermal gradient to 10°C/km. Increase the number of temperature calculations to plot from 1 to 5 and run the thermal model. Save this plot for your write-up. Is it helpful to see the temperature calculations at different times?

For each of the points above, include the requested plots with a short figure caption in your write-up, along with answers to any of the listed questions. You may want to reference these plots in your final report on these exercises, so be sure to keep the copies of the plot files.
