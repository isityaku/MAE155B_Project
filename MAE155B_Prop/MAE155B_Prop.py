### MAE155B Propulsion
# Made By: Lance De La Cruz

## Importing Libraries
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import prettytable as pt

os.system('cls') # Clearing terminal

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
prop = "APC 9x6-E"
prop_data = pd.read_csv("MAE155B_Prop/9x6E_data_9000rpm.csv")

# Configuration Parameters
B = 2 # number of blades on propeller
D = 0.2286 # outside diameter (from APC website) [in m]
S = math.pi*(D/2)**2 # disk area [in m^2]

## Propulsion Parameters
# Cruise Conditions (from APC 10x7-E Simulation Document)
V_inf = 18 # assumed cruise/freestream velocity [in m/s]
n_m = 9000 # rotational speed at cruise conditions [in RPM]
n_s = n_m/60 # rotational speed at cruise conditions [in RPS]
J_cruise = V_inf/(n_s*D) # advance ratio at cruise conditions
T_cruise = 4.825 # interpolated thrust at cruise conditions [in N]
Q_cruise = 0.130 # interpolated torque at cruise conditions [in N-m]
PWR_cruise = 122.313 # interpolated power at cruise conditions [in W]
CT_cruise = T_cruise/(rho*(n_s**2)*(D**4)) # thrust coefficient at cruise conditions
CQ_cruise = Q_cruise/(rho*(n_s**2)*(D**5)) # torque coefficient at cruise conditions
CP_cruise = PWR_cruise/(rho*(n_s**3)*(D**5)) # power coefficient at cruise conditions
eta_cruise = CT_cruise*J_cruise/CP_cruise # efficiency at cruise conditions

print("Propulsion Parameters (At Cruise Condition):")
prop_tab = pt.PrettyTable(["Parameter", "Unit", "Value"])
prop_tab.add_row(["Motor", "---", motor])
prop_tab.add_row(["Propeller", "---", prop])
prop_tab.add_row(["Number of Blades", "---", B])
prop_tab.add_row(["Outside Diameter", "m", D])
prop_tab.add_row(["Freestream Velocity", "m/s", V_inf])
prop_tab.add_row(["Propeller Rotational Velocity", "RPM", n_m])
prop_tab.add_row(["Advance Ratio", "---", J_cruise])
prop_tab.add_row(["Thrust", "N", T_cruise])
prop_tab.add_row(["Torque", "N-m", Q_cruise])
prop_tab.add_row(["Power", "W", PWR_cruise])
prop_tab.add_row(["Thrust Coefficient", "---", CT_cruise])
prop_tab.add_row(["Torque Coefficient", "---", CQ_cruise])
prop_tab.add_row(["Power Coefficient", "---", CP_cruise])
prop_tab.add_row(["Efficiency", "---", eta_cruise])

print(prop_tab)

# At variable freestream velocity
J = prop_data['J']
CT = prop_data['Thrust.1']/(rho*(n_s**2)*(D**4)) # thrust coefficient
CQ = prop_data['Torque.1']/(rho*(n_s**2)*(D**5)) # torque coefficient
CP = prop_data['PWR.1']/(rho*(n_s**3)*(D**5)) # power coefficient
eta = CT*J/CP # efficiency

## Plotting
# Thrust Coefficient v. Advance Ratio
plt.plot(J,CT)
plt.grid() # adds grid to plot
plt.xlabel("Advance Ratio (J)") # labels x axis
plt.ylabel("Thrust Coefficient") # labels y axis
plt.title("Thrust Coefficient v. Advance Ratio") # titles the plot
plt.show()

# Torque Coefficient v. Advance Ratio
plt.plot(J,CQ)
plt.grid() # adds grid to plot
plt.xlabel("Advance Ratio (J)") # labels x axis
plt.ylabel("Torque Coefficient") # labels y axis
plt.title("Torque Coefficient v. Advance Ratio") # titles the plot
plt.show()

# Power Coefficient v. Advance Ratio
plt.plot(J,CP)
plt.grid() # adds grid to plot
plt.xlabel("Advance Ratio (J)") # labels x axis
plt.ylabel("Power Coefficient") # labels y axis
plt.title("Power Coefficient v. Advance Ratio") # titles the plot
plt.show()

# Efficiency v. Advance Ratio
plt.plot(J,eta)
plt.grid() # adds grid to plot
plt.xlabel("Advance Ratio (J)") # labels x axis
plt.ylabel("Efficiency") # labels y axis
plt.title("Efficiency v. Advance Ratio") # titles the plot
plt.show()
