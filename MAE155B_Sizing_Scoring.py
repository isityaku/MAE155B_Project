### MAE155B Sizing and Scoring Analysis
# Made By: Lance De La Cruz

## Outline
# Define the following parameters:
#   Atmospheric Parameters
#   Propulsion Parameters
#   Sizing / Weight Parameters
#   Aerodynamic Parameters
# Conduct Preliminary Analysis for the following:
#   Stall
#   Climb
#   Cruise
#   Maneuver
#   Takeoff
#   Landing
# Create the sizing plot
# Use sizing plot to determine flight score by obtaining the following:
#   Wing Loading
#   Thrust-to-Weight Ratio
#   Gross and Payload Ratio
#   Wing Area, Span, and Chord
#   Planform Dimension Sum
#   Flight Score

## Importing Libraries
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Clearing terminal
os.system('cls')

## Defining Parameters
# Atmospheric Parameters
g = 9.81 # gravitational acceleration [in m/s^2]
gamma = 1.4 # specific heat ratio
P = 101280 # atmospheric pressure [in Pa]
temp = 286; # atmospheric temperature [in K]
R = 287 # gas constant
mu = 1.8*10**-5 # dynamic viscosity
rho = P/(R*temp) # atmospheric density [in kg/m^3]
a = math.sqrt(gamma*R*temp) # speed of sound [in m/s]

# Propulsion Parameters
motor = "Cobra C-2217/16 Brushless Motor"
prop = "APC 11x4.7-SF"
T = 0.879*9.81 # motor/prop configuration max thrust [in N]

# Sizing Parameters
W_S = 9.76*9.81 # assumed wing loading [in N/m^2]
T_W = 0.8 # assumed thrust-to-weight ratio
AR = 6 # assumed aspect ratio
B = 1 # empty weight bonus
p_weight = 0.037*9.81 # individual payload weight [in N/syringe]
Wp = p_weight*np.arange(0,50,1) # payload weight (ranging from 0 to 100 units of payload) [in N]
Wp_W0 = 0.4 # assuemd payload-gross weight fraction
L_b = 1/1.5 # length-to-width ratio (length-to-wingspan)

# Aerodynamic Parameters
CL_max = 1 # assumed max lift coefficient
CL_cruise = 0.5 # assumed cruise lift coefficient
Cd0 = 0.04 # assumed zero AoA drag coefficient
e = 1.78*(1-(0.045*AR**0.68))-0.64 # oswald efficiency factor
n = 3.5 # max load factor (used for maneuver calculation)

## Preliminary Analysis
# Stall
V_stall = math.sqrt(2*W_S/(rho*CL_max)) # Stall Velocity

# Climb
gamma = 15 # assumed climb angle [in degrees]
V_climb = 8 # climb velocity
G = math.sin(math.radians(gamma)) # climb gradient
q_climb = 0.5*rho*V_climb**2 # climb dynamic pressure

# Cruise
V_cruise = 18 # assumed cruise velocity [in m/s]
q_cruise = 0.5*rho*V_cruise**2 # cruise dynamic pressure

# Ceiling
D_L_min = 2*math.sqrt(Cd0/(math.pi*e*AR)) # minimum drag-lift ratio

## Sizing Plot
# Sizing Equations
W_S_plot = np.linspace(0,500,1000)
T_W_climb = (G+(W_S_plot/(math.pi*e*AR*q_climb))+(Cd0*q_climb/W_S_plot)) # climb equation
T_W_ceil = (0*W_S_plot)+D_L_min # ceiling equation
T_W_man = (Cd0*q_cruise/W_S_plot) + ((n**2)*W_S_plot/(math.pi*e*AR*q_cruise)) # maneuver equation

# Plotting Functions
plt.plot(W_S_plot, T_W_climb, label="Climb", color="blue") # plots climb equation
plt.plot(W_S_plot, T_W_ceil, label="Ceiling", color="red") # plots ceiling equation
plt.plot(W_S_plot, T_W_man, label="Maneuver", color="green") # plots maneuver equation
plt.vlines(W_S, 0, 1, label="Stall", color="black") # plots stall wing loading
plt.xlim(0, 150) # adjusts x axis
plt.ylim(0, 1) # adjusts y axis
plt.grid() # adds grid to plot
plt.xlabel("W/S [in N/m^2]") # labels x axis
plt.ylabel("T/W") # labels y axis
plt.title("Thrust-to-Weight v. Wing Loading [in N/m^2]") # titles the plot
plt.legend() # displays the legend
plt.show() # produces figure

## Flight Score Calculation
W_S_new = 82.4756 # wing loading from sizing plot
T_W_new = 0.405472 # thrust-to-weight ratio from sizing plot
W0_frac = Wp*(1/Wp_W0) # gross weight based on payload-gross weight fraction [in N]
W0_T = T*(1/T_W_new) # gross weight based on max thrust from motor-prop config. [in N]
W0_temp = np.where(W0_frac < W0_T)[0] # finds array positions for gross weight values from weight fraction less than gross weight based on T/W
W0 = W0_frac[np.max(W0_temp)] # gross weight [in N]
Wp = W0*Wp_W0 # payload weight [in N]
p = Wp/p_weight # units of payload [in syringes]
S = W0*(1/W_S_new) # wing area [in m^2]
b = math.sqrt(AR*S) # wing span [in m]
c = S/b # wing chord [in m]
D = math.sqrt(AR*W0/W_S_new)*(1+L_b) # sum of aircraft planform dimensions [in m]
FS = 20*Wp + Wp_W0 - D**4 + 20*B # flight score function