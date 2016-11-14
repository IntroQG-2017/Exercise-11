# river_profiles.py
#
# This script calculates bedrock river profiles as a function of rock uplift
# rate and various stream-power erosion law parameters.
#
# User-defined variables are listed at the top of the script.
#
# This code is based on a similar MATLAB code written by Brian Yanites and
# Todd Ehlers.
#
# Dave Whipp - 13.04.16

# Import NumPy, Matplotlib and system libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys

# Define function to update river profile elevations
def update_topo(topography, dhdt, dx, K, area, m, n, uplift, dt, dtplot):
    for step in range(int(dtplot/dt)):           # Loop over amount of time between plots
        slope = abs(np.diff(topography))/dx      # Calculate profile slopes; diff() calculates elevation difference between points
        slope = np.append(slope, slope[-1])      # Append one extra slope value to be same size as area array
        dhdt = -K * area**m * slope**n           # Calculate stream-power erosion rate
        topography += (dhdt + uplift) * dt       # Change topography based on uplift and river erosion
        tfix = np.where(topography < 0.0)        # Find any values of topography less than 0.0
        topography[tfix] = 0.0                   # Set those values equal to 0.0
    return topography, dhdt                      # Return updated elevations and erosion rates for river profile

# Define function to animate river profile elevations
def animate(i, topography, dhdt, dx, K, area, m, n, uplift, dt, dtplot, plottime):
    topography,dhdt=update_topo(topography, dhdt, dx, K, area, m, n, uplift,    # Update elevations on river profile
                               dt, dtplot)
    timenow = plottime[i]                                      # Set timenow to current model time
    elevmax = max(topography)                                  # Calculate maximum elevation of river profile
    plot1.set_ydata(topography)                                # Update plot topography
    axis1.set_ylim([0.0, max(topography) * 1.1])               # Update y-axis range for new topography
    timetext.set_y(max(topography) * 0.9)                      # Reposition time text
    timetext.set_text(str(timenow)+" years")                   # Update time text
    elevtext.set_y(max(topography) * 0.8)                      # Reposition time text
    elevtext.set_text("Max elev: {0:.1f} m".format(elevmax))   # Update time text
    # ADD STUFF TO UPDATE SECOND SUBPLOT HERE
    #
    #
    return plot1,                                              # Return updated plot

# Main program code below

#--- USER-DEFINED VARIABLES BELOW -----------------------------------------#
max_elevation = 1500.0      # Maximum elevation for initial topography [m]
min_elevation = 0.0         # Minimum elevation for initial topography [m]
x_length = 100000.0         # Lateral length of river profile [m]
dx = 1000.0                 # Spacing of points along profile for elevation calculations [m]
topo_option = 2             # Flag for initial topography; 1 = Constant slope, 2 = Flat
simulation_time = 2000000.0 # Total simulation time [a]
dt = 2.0                    # Time step size for profile calculations [a]
dtplot = 2000.0            # Time step for displaying plots [a]

# EROSIONAL CONSTANTS
m = 1.0               # Area exponent
                      #   m = 0.3 for bed shear stress representation
                      #   m = 1.0 for stream power per unit channel length
                      #   m = 0.5 for stream power per unit bed area
n = 2.0               # Slope exponent
                      #   n = 0.7 for bed shear stress representation
                      #   n = 1.0 for stream power per unit channel length
                      #   n = 1.0 for stream power per unit bed area
h = 1.69              # Scaling exponent for distance-to-area conversion
ka = 6.69             # Distance/area scaling coefficient for distance-to-area conversion [m**0.33]
K = 1.6E-8            # Fluvial erosional efficiency factor [m**0.33 / a]
                      # This factor accounts for lithology, climate, channel geometry, and perhaps sediment supply (!)
uplift_rate = 0.001   # Rock uplift rate [m/a]

#--- END USER-DEFINED VARIABLES -------------------------------------------#

# Define arrays needed for river profile calculations
xpoints = x_length/dx                         # Number of points for x array
x = np.linspace(0.0, x_length, xpoints + 1)   # Define x array from 0 to L by dx
uplift = np.zeros(len(x)) + uplift_rate       # Create rock uplift array with constant uplift rate
xkm = x / 1000.0                              # Create xkm array with x values converted to km
dhdt = np.zeros(len(x))                       # Erosion rate array, same size as x

# SELECT INPUT TOPOGRAPHY
# OPTION 1: Initial topography is a sloping surface
if topo_option == 1:
    topography = np.linspace(max_elevation, min_elevation, xpoints + 1)
# OPTION 2: Initial topography is a flat surface
elif topo_option == 2:
    topography = np.zeros(len(x)) + max_elevation
    topography[-1] = min_elevation
# Stop with error message if given a bad value
else:
    print("Bad value listed for topo_option. topo_option must be '1' or '2'.")
    sys.exit()

# Calculate drainage basin area as a function of distance from divide x
area = ka*x**h
area[0] = area[1]

# Fill values for time and plottime arrays
time = np.linspace(0.0, simulation_time, simulation_time / dt + 1)
plottime = np.linspace(0.0, simulation_time, simulation_time / dtplot + 1)

# Create plotting window
fig = plt.figure()

# Format subplot 1
axis1 = plt.subplot(1,1,1)                  # Set axis1 as the first plot
axis1.set_xlim([0.0, max(xkm)])             # Set the x-axis limits for plot 1
axis1.set_ylim([0.0, max_elevation*1.1])    # Set the y-axis limits for plot 1
plot1, = plt.plot(xkm, topography)          # Define plot1 as the first plot
plt.xlabel("Distance from drainage divide [km]")
plt.ylabel("River channel elevation [m]")
plt.title("River channel profile evolution model")
# Define time for display on plot
timenow = 0.0
timetext = plt.text(70.,max_elevation*0.9, str(timenow)+" years")
elevmax = max(topography)
elevtext = plt.text(70., max_elevation*0.8, "Max elev: {0:.1f} m".format(elevmax))

# Format subplot 2
#
# FIGURE 2 STUFF GOES HERE
#

# Animate and display plot results
anim = animation.FuncAnimation(fig, animate, range(len(plottime)),
                              fargs=(topography, dhdt, dx, K, area, m, n,
                              uplift, dt, dtplot, plottime), interval=20,
                              blit=False, repeat=False)
plt.show()
