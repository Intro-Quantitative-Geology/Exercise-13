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
plotT0 = True          # Plot the initial thermal solution
plotTzHistory = False  # Plot the temperature-depth history of the tracked particle
calcAHe = False        # Calculate apatite (U-Th)/He age
calcZHe = False        # Calculate zircon (U-Th)/He age
calcMAr = False        # Calculate muscovite 40Ar/39Ar age

# THERMAL MODEL PARAMETERS
Tgradient = 10.0  # Initial thermal gradient [deg. C / km]
zMax = 50.0       # Thermal model thickness [km]
tMax = 50.0       # Total thermal model simulation time [Ma]
vz = 0.5          # Vertical advection velocity [km/Ma]
kappa = 32.0      # Thermal diffusivity [km2 / Ma]
npz = 101         # Number of depth points for temperature calculation
npt = 401         # Number of times for temperature calculation
highT = 1000.0    # Temperature to assign if tracked particle depth exceeds zmax
nTout = 1         # Number of temperature profiles to plot

# CLOSURE TEMPERATURE PARAMETERS
# Apatite (U-Th)/He
aAHe = 100.0       # Apatite grain size [um]
EaAHe = 138.0e3    # Activation energy [J/mol]
AAHe = 25.0        # Geometry factor [25 for sphere]
D0AHe = 5.0e-3     # Diffusivity at infinite temperature [m2/s]
tauAHe = 1.0       # Initial guess for characteristic time

# Zircon (U-Th)/He
aZHe = 100.0       # Zircon grain size [um]
EaZHe = 168.0e3    # Activation energy [J/mol]
AZHe = 25.0        # Geometry factor [25 for sphere]
D0ZHe = 4.6e-5     # Diffusivity at infinite temperature [m2/s]
tauZHe = 1.0       # Initial guess for characteristic time

# Muscovite Ar/Ar
aMAr = 500.0       # Muscovite grain size [um]
EaMAr = 183.0e3    # Activation energy [J/mol]
AMAr = 8.7         # Geometry factor [8.7 for planar sheet]
D0MAr = 3.3e-6     # Diffusivity at infinite temperature [m2/s]
tauMAr = 1.0       # Initial guess for characteristic time

# OTHER CONSTANTS
R = 8.314         # Universal gas constant

#--- END USER-DEFINED VARIABLES -----------------------------------------------#

# Set initial thermochronometer ages
ageA = tMax
ageZ = tMax
ageM = tMax

# Convert units
aAHe = aAHe / 1.0e6 / 1.0e3                                            # um -> km
aZHe = aZHe / 1.0e6 / 1.0e3                                            # um -> km
aMAr = aMAr / 1.0e6 / 1.0e3                                            # um -> km
D0AHe = D0AHe * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma
D0ZHe = D0ZHe * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma
D0MAr = D0MAr * (1 / 1000.0**2.0) * (1.0e6 * 365.25 * 24.0 * 3600.0)   # m2/s -> km2/Ma

# Calculate diffusion parameter D0/a2
D0a2AHe = D0AHe / aAHe**2.0    # Apatite
D0a2ZHe = D0ZHe / aZHe**2.0    # Zircon
D0a2MAr = D0MAr / aMAr**2.0    # Muscovite

# Thermal model setup
z = np.linspace(0.0, zMax, npz)            # Define depth range array
t = np.linspace(0.0, tMax, npt)            # Define time range array
tMa = np.linspace(tMax, 0.0, npt)          # Define time array in Ma (time before present)
zHistory = np.linspace(tMax * vz, 0.0, npt)   # Define z particle position history
Thistory = np.zeros(len(t))                   # Define initial temperature history array
Temperature = np.zeros(len(z))                # Define initial temperature array

iout = int(float(len(t)-1)/float(nTout))   # Define increment for ploting output

# Make a new plot window
plt.figure()

# Loop over all times and calculate temperature
for i in range(len(t)):
    # Loop over all depths and calculate temperature
    for j in range(len(z)):
        # Call function to calculate temperature for this depth and time
        Temperature[j] = transientTemp1D(z[j],t[i],Tgradient,vz,kappa)
    # Set the temperature history temperature to highT if the tracked particle
    # depth is below zmax (where temperature would be undefined)
    if zHistory[i] > max(z):
        Thistory[i] = highT
    # Otherwise, store the current temperature at the depth of the tracked particle
    else:
        Thistory[i] = transientTemp1D(zHistory[i],t[i],Tgradient,vz,kappa)

    # If plotting of the initial geotherm is requested, make the plot
    if plotT0 and i == 0:
        plt.plot(Temperature,-z,label=str(tMa[i])+" Ma")

    # If the current iteration is one of the plotting increments, make the plot
    if i == iout:
        iout = iout + int(float(len(t))/float(nTout))
        plt.plot(Temperature,-z,label=str(tMa[i])+" Ma")

# Set the initial closure temperatures for the previous step to one
TcAHeP = 1.0
TcZHeP = 1.0
TcMArP = 1.0

# Set the previous temperature to the max value in the current temperature array
ThistoryP = max(Temperature)

# Loop over all positions in the temperature history array
for i in range(len(Thistory)):
    # Calculate the cooling rate for the first temperature value
    if i == 0:
        dTdt = (Thistory[i] - Thistory[i+1]) / (tMa[i] - tMa[i+1])
    # Calculate the cooling rate for the last temperature value
    elif i == len(Thistory)-1:
        dTdt = (Thistory[i-1] - Thistory[i]) / (tMa[i-1] - tMa[i])
    # Calculate the cooling rate for the the intermediate temperature values
    else:
        dTdt = (Thistory[i-1] - Thistory[i+1]) / (tMa[i-1] - tMa[i+1])

    # Ensure the cooling rate is at least 1 deg. C per 10 Ma
    dTdt = max(dTdt, 0.1 / (1.0e6 * 365.25 * 24.0 * 3600.0))

    # Calculate apatite (U-Th)/He closure temperature if requested
    if calcAHe:
        # Calculate diffusivity characteristic time
        tauAHe = (R * (Thistory[i] + 273.15)**2.0) / (EaAHe * dTdt)
        # Calculate new closure temperature
        TcAHe = EaAHe / (R * np.log(AAHe * tauAHe * D0a2AHe)) - 273.15
        # Calculate new cooling age
        if Thistory[i] > TcAHe:
            ratio = (TcAHeP - ThistoryP)/(TcAHeP - ThistoryP + Thistory[i] - TcAHe)
            ageA = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        TcAHeP = TcAHe

    # Calculate zircon (U-Th)/He closure temperature if requested
    if calcZHe:
        # Calculate diffusivity characteristic time
        tauZHe = (R * (Thistory[i] + 273.15)**2.0) / (EaZHe * dTdt)
        # Calculate new closure temperature
        TcZHe = EaZHe / (R * np.log(AZHe * tauZHe * D0a2ZHe)) - 273.15
        # Calculate new cooling age
        if Thistory[i] > TcZHe:
            ratio = (TcZHeP - ThistoryP)/(TcZHeP - ThistoryP + Thistory[i] - TcZHe)
            ageZ = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        TcZHeP = TcZHe

    # Calculate muscovite Ar/Ar closure temperature if requested
    if calcMAr:
        # Calculate diffusivity characteristic time
        tauMAr = (R * (Thistory[i] + 273.15)**2.0) / (EaMAr * dTdt)
        # Calculate new closure temperature
        TcMAr = EaMAr / (R * np.log(AMAr * tauMAr * D0a2MAr)) - 273.15
        # Calculate new cooling age
        if Thistory[i] > TcMAr:
            ratio = (TcMArP - ThistoryP)/(TcMArP - ThistoryP + Thistory[i] - TcMAr)
            ageM = tMa[i] + (tMa[i-1] - tMa[i]) * ratio
        TcMArP = TcMAr

    # Store previous temperature in thermal history
    ThistoryP = Thistory[i]

# Write apatite (U-Th)/He age to screen if requested
if calcAHe:
    print("Apatite (U-Th)/He age: "+str(ageA))

# Write zircon (U-Th)/He age to screen if requested
if calcZHe:
    print("Zircon (U-Th)/He age: "+str(ageZ))

# Write muscovite Ar/Ar age to screen if requested
if calcMAr:
    print("Muscovite Ar/Ar age: "+str(ageM))

# Plot particle depth-temperature history if requested
if plotTzHistory:
    plt.plot(Thistory,-zHistory,'*')

# Label axes and add title
plt.xlabel("")
plt.ylabel("")
plt.title("")

# Display line legend
plt.legend()

# Display figure
plt.show()
