"""
age_predict_1D.py

This code predicts thermochronometer cooling ages using a 1-D transient
advection-diffusion thermal model and Dodson's method for predicting mineral
closure temperatures.

The surface temperature is assumed to be 0 deg. C and the final simulation
time in the model is 0 Ma.

@author: Dave Whipp - 5.12.16
"""

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# Define function to calculate temperature as a function of depth, time, initial
# thermal gradient, advection velocity and thermal diffusivity
def transientTemp1D(z,t,G,vz,kappa):
    # Calculate T separately for case where t = 0 to avoid divide-by-zero warnings
    if t == 0:
        Temperature = G * z
    else:
        Temperature = G * (z + vz * t) + (G / 2.0) * ((z - vz * t) * np.exp(-(vz * z) / kappa) *\
        erfc((z - vz * t) / (2.0 * np.sqrt(kappa * t))) - (z + vz * t) * \
        erfc((z + vz * t) / (2.0 * np.sqrt(kappa * t))))
    return Temperature

#--- USER-DEFINED VARIABLES ---------------------------------------------------#

# PROGRAM OPTIONS
plot_T0 = True       # Plot the initial thermal solution
plot_Tzhist = False  # Plot the temperature-depth history of the tracked particle
calc_AHe = False     # Calculate apatite (U-Th)/He age
calc_ZHe = False     # Calculate zircon (U-Th)/He age
calc_MAr = False     # Calculate muscovite 40Ar/39Ar age

# THERMAL MODEL PARAMETERS
Tgrad = 10.0      # Initial thermal gradient [deg. C / km]
zmax = 50.0       # Thermal model thickness [km]
tmax = 50.0       # Total thermal model simulation time [Ma]
vz = 0.5          # Vertical advection velocity [km/Ma]
kappa = 32.0      # Thermal diffusivity [km2 / Ma]
npz = 101         # Number of depth points for temperature calculation
npt = 401         # Number of times for temperature calculation
highT = 1000.0    # Temperature to assign if tracked particle depth exceeds zmax
nTout = 1         # Number of temperature profiles to plot

# CLOSURE TEMPERATURE PARAMETERS
# Apatite (U-Th)/He
a_a = 100.0       # Apatite grain size [um]
Ea_a = 138.0e3    # Activation energy [J/mol]
A_a = 25.0        # Geometry factor [25 for sphere]
D0_a = 5.0e-3     # Diffusivity at infinite temperature [m2/s]
tau_a = 1.0       # Initial guess for characteristic time

# Zircon (U-Th)/He
a_z = 100.0       # Zircon grain size [um]
Ea_z = 168.0e3    # Activation energy [J/mol]
A_z = 25.0        # Geometry factor [25 for sphere]
D0_z = 4.6e-5     # Diffusivity at infinite temperature [m2/s]
tau_z = 1.0       # Initial guess for characteristic time

# Muscovite Ar/Ar
a_m = 500.0       # Muscovite grain size [um]
Ea_m = 183.0e3    # Activation energy [J/mol]
A_m = 8.7         # Geometry factor [8.7 for planar sheet]
D0_m = 3.3e-6     # Diffusivity at infinite temperature [m2/s]
tau_m = 1.0       # Initial guess for characteristic time

# OTHER CONSTANTS
R = 8.314         # Universal gas constant

#--- END USER-DEFINED VARIABLES -----------------------------------------------#

# Set initial thermochronometer ages
ageA = tmax
ageZ = tmax
ageM = tmax

# Convert units
a_a = a_a / 1.0e6 / 1.0e3                                            # um -> km
a_z = a_z / 1.0e6 / 1.0e3                                            # um -> km
a_m = a_m / 1.0e6 / 1.0e3                                            # um -> km
D0_a = D0_a * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma
D0_z = D0_z * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma
D0_m = D0_m * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma

# Calculate diffusion parameter D0/a2
D0a2_a = D0_a / a_a**2.0    # Apatite
D0a2_z = D0_z / a_z**2.0    # Zircon
D0a2_m = D0_m / a_m**2.0    # Muscovite

# Thermal model setup
z = np.linspace(0.0, zmax, npz)            # Define depth range array
t = np.linspace(0.0, tmax, npt)            # Define time range array
tMa = np.linspace(tmax, 0.0, npt)          # Define time array in Ma (time before present)
zhist = np.linspace(tmax * vz, 0.0, npt)   # Define z particle position history
Thist = np.zeros(len(t))                   # Define initial temperature history array
T = np.zeros(len(z))                       # Define initial temperature array

iout = int(float(len(t)-1)/float(nTout))   # Define increment for ploting output

# Make a new plot window
plt.figure()

# Loop over all times and calculate temperature
for i in range(len(t)):
    # Loop over all depths and calculate temperature
    for j in range(len(z)):
        # Call function to calculate temperature for this depth and time
        T[j] = transient_temp_1D(z[j],t[i],Tgrad,vz,kappa)
    # Set the temperature history temperature to highT if the tracked particle
    # depth is below zmax (where temperature would be undefined)
    if zhist[i] > max(z):
        Thist[i] = highT
    # Otherwise, store the current temperature at the depth of the tracked particle
    else:
        Thist[i] = transient_temp_1D(zhist[i],t[i],Tgrad,vz,kappa)

    # If plotting of the initial geotherm is requested, make the plot
    if plot_T0 and i == 0:
        plt.plot(T,-z,label=str(tMa[i])+" Ma")

    # If the current iteration is one of the plotting increments, make the plot
    if i == iout:
        iout = iout + int(float(len(t))/float(nTout))
        plt.plot(T,-z,label=str(tMa[i])+" Ma")

# Set the initial closure temperatures for the previous step to one
Tc_ap = 1.0
Tc_zp = 1.0
Tc_mp = 1.0

# Set the previous temperature to the max value in the current temperature array
Thistp = max(T)

# Loop over all positions in the temperature history array
for i in range(len(Thist)):
    # Calculate the cooling rate for the first temperature value
    if i == 0:
        dTdt = (Thist[i] - Thist[i+1]) / (tMa[i] - tMa[i+1])
    # Calculate the cooling rate for the last temperature value
    elif i == len(Thist)-1:
        dTdt = (Thist[i-1] - Thist[i]) / (tMa[i-1] - tMa[i])
    # Calculate the cooling rate for the the intermediate temperature values
    else:
        dTdt = (Thist[i-1] - Thist[i+1]) / (tMa[i-1] - tMa[i+1])

    # Ensure the cooling rate is at least 1 deg. C per 10 Ma
    dTdt = max(dTdt, 0.1 / (1.0e6 * 365.25 * 24.0 * 3600.0))

    # Calculate apatite (U-Th)/He closure temperature if requested
    if calc_AHe:
        # Calculate diffusivity characteristic time
        tau_a = (R * (Thist[i] + 273.15)**2.0) / (Ea_a * dTdt)
        # Calculate new closure temperature
        Tc_a = Ea_a / (R * np.log(A_a * tau_a * D0a2_a)) - 273.15
        # Calculate new cooling age
        if Thist[i] > Tc_a:
            ratio = (Tc_ap - Thistp)/(Tc_ap - Thistp + Thist[i] - Tc_a)
            ageA = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        Tc_ap = Tc_a

    # Calculate zircon (U-Th)/He closure temperature if requested
    if calc_ZHe:
        # Calculate diffusivity characteristic time
        tau_z = (R * (Thist[i] + 273.15)**2.0) / (Ea_z * dTdt)
        # Calculate new closure temperature
        Tc_z = Ea_z / (R * np.log(A_z * tau_z * D0a2_z)) - 273.15
        # Calculate new cooling age
        if Thist[i] > Tc_z:
            ratio = (Tc_zp - Thistp)/(Tc_zp - Thistp + Thist[i] - Tc_z)
            ageZ = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        Tc_zp = Tc_z

    # Calculate muscovite Ar/Ar closure temperature if requested
    if calc_MAr:
        # Calculate diffusivity characteristic time
        tau_m = (R * (Thist[i] + 273.15)**2.0) / (Ea_m * dTdt)
        # Calculate new closure temperature
        Tc_m = Ea_m / (R * np.log(A_m * tau_m * D0a2_m)) - 273.15
        # Calculate new cooling age
        if Thist[i] > Tc_m:
            ratio = (Tc_mp - Thistp)/(Tc_mp - Thistp + Thist[i] - Tc_m)
            ageM = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        Tc_mp = Tc_m

    # Store previous temperature in thermal history
    Thistp = Thist[i]

# Write apatite (U-Th)/He age to screen if requested
if calc_AHe:
    print("Apatite (U-Th)/He age: "+str(ageA))

# Write zircon (U-Th)/He age to screen if requested
if calc_ZHe:
    print("Zircon (U-Th)/He age: "+str(ageZ))

# Write muscovite Ar/Ar age to screen if requested
if calc_MAr:
    print("Muscovite Ar/Ar age: "+str(ageM))

# Plot particle depth-temperature history if requested
if plot_Tzhist:
    plt.plot(Thist,-zhist,'*')

# Label axes and add title
plt.xlabel("")
plt.ylabel("")
plt.title("")

# Display line legend
plt.legend()

# Display figure
plt.show()
