### MAE155B Propulsion [WIP]
# Made By: Lance De La Cruz

## Outline
# State and justify motor-prop configuration.
# Find the following parameters (at cruise):
#   Efficiency
#   Thrust
#   Torque
#   Power
#   Advance Ratio
#   Rotational Speed
# Plot the following v. advance ratio:
#   Efficiency
#   Thrust Coefficient
#   Torque Coefficient
#   Power Coefficient

## Importing Libraries
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Atmospheric Parameters
g = 9.81 # gravitational acceleration [in m/s^2]
gamma = 1.4 # specific heat ratio
P = 101280 # atmospheric pressure [in Pa]
temp = 286; # atmospheric temperature [in K]
R = 287 # gas constant
mu = 1.8*10**-5 # dynamic viscosity
rho = P/(R*temp) # atmospheric density [in kg/m^3]
a = math.sqrt(gamma*R*temp) # speed of sound [in m/s]

## Propulsion Configuration (refer to MAE155B_Sizing_Scoring.py)
# Motor/Prop Configuration
motor = "Cobra C-2217/16 Brushless Motor"
prop = "APC 10x7-E"

# Configuration Parameters
T = 0.687*9.81 # motor/prop configuration max thrust [in N]
Q = [] # motor/prop configuration torque
Power = [] # motor/prop configuration power
B = 2 # number of blades on propeller
D = 0.2794 # outside diameter (from APC website) [in m]
S = math.pi*(D/2)**2 # disk area [in m^2]
T_S = T/S # disk loading [in N/m^2]
pitch = 0.1778 # propeller pitch [in m]
beta = math.atan2(pitch/(2*math.pi*(D/2)))

## Propulsion Parameters
# At cruise conditions
V_inf = 18 # assumed cruise/freestream velocity [in m/s]
alpha = 0 # angle-of-attack (theoretically zero for cruise)
phi = beta-alpha
w_cruise = V_inf/((D/2)*math.tan(math.radians(phi))) # angular speed at cruise [rad/s]
n_cruise = w_cruise*0.159 # rotational speed at cruise [rev/s]
J_cruise = V_inf/(n_cruise*D) # advance ratio at cruise
nu_cruise = T*V_inf/(Q*w_cruise) # efficiency at cruise
CT_cruise = T/(rho*(n_cruise**2)*(D**4)) # thrust coefficient at cruise
CQ_cruise = Q/(rho*(n_cruise**2)*(D**5)) # torque coefficient at cruise
CP_cruise = Power/(rho*(n_cruise**3)*(D**5)) # power coefficient at cruise

# At variable freestream velocity
V = np.linspace(0,500,500) # variable freestream velocity
w = V/((D/2)*math.tan(math.radians(phi))) # angular speed [rad/s]
n = w*0.159 # rotational speed
J = V_inf/(n*D) # advance ratio
CT = T/(rho*(n**2)*(D**4)) # thrust coefficient
CQ = Q/(rho*(n**2)*(D**5)) # torque coefficient
CP = Power/(rho*(n**3)*(D**5)) # power coefficient
nu = CT*J/CP # efficiency