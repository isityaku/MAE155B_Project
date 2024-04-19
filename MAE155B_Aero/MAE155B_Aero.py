### MAE155B Aerodynamics [WIP]
# Made By: Lance De La Cruz

## Outline
# State and justify airfoil choice
#   Include lift and drag polars for airfoil (XFLR5)
# Compute parameters for aircraft lift and drag coefficient equations
#   Refer to lecture for parameters required
# Define the following aerodynamic parameters (at cruise):
#   Lift Coefficient
#   Drag Coefficient
#   Lift-to-Drag Ratio

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
L = 0.75 # estimated aircraft planform length (refer to MAE155B_Sizing_Scoring.py) [in m]

## Cruise Conditions
V_cruise = 18 # assumed cruise velocity [in m/s]
M = V_cruise/a # cruise mach number
Re = rho*V_cruise*L/mu # Reynold's number at cruise velocity

## Airfoil Configuration
# Wing
foil_wing = "CLARK Y" # airfoil used for aircaft wing
foil_wing_data = pd.read_csv("MAE155B_Data/xf-clarky-il-1000000.csv")
Cl_alp = np.polyfit(foil_wing_data["Alpha"], foil_wing_data["Cl"], 1)[0]
alp0_Cl = foil_wing_data.loc[foil_wing_data["Cl"] == foil_wing_data["Cl"].abs().min()] # lift coefficient = 0
alp0 = alp0_Cl["Alpha"].max() # angle-of-attack where lift coefficient = 0
Cl_max = foil_wing_data["Cl"].max() # maximum 2-D lift coefficient
Cl_min = foil_wing_data["Cl"].min() # minimum 2-D lift coefficient
Cd_max = foil_wing_data["Cd"].max() # maximum 2-D drag coefficient
Cd_min = foil_wing_data["Cd"].min() # minimum 2-D drag coefficient

# Horizontal + Vertical Stabilizer
foil_tail = "NACA 0015" # airfoil used for aircraft tail

## Aircraft Configuration (refer to MAE155B_Sizing_Scoring.py)
# Wing
b_wing = 1.359 # wing span [in m]
c_wing = 0.227 # wing chord [in m]
t_c_wing = 0.117 # wing thickness-to-chord ratio
x_c_wing = 0.28 # location of maximum thickness for wing
MAC = (2/3)*c_wing # mean aerodynamic chord [in m]
S_ref_wing = b_wing*c_wing # wing area [in m^2]
AR = b_wing/c_wing # aspect ratio
e = 1.78*(1-(0.045*AR**0.68))-0.64 # oswald efficiency factor
d = 0.2 # fuselage diameter [in m]
S_exp_ref = 0.85 # exposed-reference planform area ratio
S_wet_wing = 0.95*S_ref_wing # wetted area for wing
Re_wing = rho*V_cruise*b_wing/mu
Cf_wing = 0.455/((math.log(Re_wing, 10)**2.58)*(1+(0.144*M**2))**0.65) # turbulent skin friction coefficient

# Horizontal Stabilizer
b_hs = 0.5 # assumed horizontal stabilizer span [in m]
c_hs = 0.125 # assumed horizontal stabilizer chord lenght [in m]
S_ref_hs = b_hs*c_hs # horizontal stabilizer area
t_c_hs = 0.15 # horizontal stabilizer thickness-to-chord ratio
x_c_hs = 0.28 # location of maximum thickness for horizontal stabilizer
S_wet_hs = 0.95*S_ref_hs # wetted area for horizontal stabilizer
Re_hs = rho*V_cruise*b_hs/mu
Cf_hs = 0.455/((math.log(Re_hs, 10)**2.58)*(1+(0.144*M**2))**0.65) # turbulent skin friction coefficient

# Vertical Stabilizer
b_vs = 0.25 # assumed vertical stabilizer span [in m]
c_vs = 0.125 # assumed vertical stabilizer chord lenght [in m]
S_ref_vs = b_vs*c_vs # vertical stabilizer area
t_c_vs = 0.15 # vertical stablizer thickness-to-chord ratio
x_c_vs = 0.28 # location of maximum thickness for vertical stabilizer
S_wet_vs= 0.95*S_ref_vs # wetted area for vertical stabilizer
Re_vs = rho*V_cruise*b_vs/mu
Cf_vs = 0.455/((math.log(Re_vs, 10)**2.58)*(1+(0.144*M**2))**0.65) # turbulent skin friction coefficient

# Fuselage
l_fuse = 0.75 # assumed fuselage length (planform) [in m]
w_fuse = 0.25 # assumed fuselage width (frontal) [in m]
h_fuse = 0.2 # assumed fuselage height (frontal) [in m]
A_fuse = w_fuse*h_fuse # frontal area of fuselage [in m2]
f_fuse = l_fuse/math.sqrt((4/math.pi)*A_fuse) # maximum thickness location for fuselage [in m]
S_wet_fuse = (2*w_fuse*h_fuse)+(2*l_fuse*w_fuse)+(2*l_fuse*h_fuse) # wetted area for fuselage
Re_fuse = rho*V_cruise*l_fuse/mu
Cf_fuse = 0.455/((math.log(Re_fuse, 10)**2.58)*(1+(0.144*M**2))**0.65) # turbulent skin friction coefficient

# Landing Gear
Cd_land = 1.01 # assumed drag coefficient for landing gear

# Motor
Cd_motor = 0.34 # assumed drag coefficient for motor
d_motor = 0.0277 # outside diameter of motor [in m]
A_motor = math.pi*(d_motor/2)**2 # frontal area of the motor [in m2]

# Misc
Lambda_quart = 0 # sweep angle at quarter-chord
Lambda_max = 0 # max sweep angle
Lambda_m = 0 # some sort of sweep angle (doesn't rlly matter, only goes to trig functions)

## 3-D Coefficient Calculation
alp = np.linspace(-5,20,100)

# 3-D Lift Calculation
beta = math.sqrt(1-(M**2))
F = 1.07*(1+(d/b_wing))**2
eta = (Cl_alp*180/math.pi)/(2*math.pi/beta)
CL_alp = (2*math.pi*AR)/(2+(math.sqrt(4+((AR**2)*(beta**2)/(eta**2))*(1+(math.tan(math.radians(Lambda_max))**2/(beta**2))))))*S_exp_ref*F*math.pi/180 # semi-empirical formula for 3-D lift-curve slope
CL = CL_alp*(alp-alp0)
CL_cruise = CL_alp*(-alp0)
CL_max = 0.9*Cl_max*math.cos(Lambda_quart)
CL_min = min(CL)

# 3-D Drag Calculation
K = np.polyfit((foil_wing_data["Cl"]-Cl_min)**2, foil_wing_data["Cd"], 1)[0]
Q = 1.1
FF_wing = (1+(0.6*t_c_wing/x_c_wing)+(100*t_c_wing**4))*((1.34*M**0.18)*(math.cos(Lambda_m)**0.28))
FF_hs = (1+(0.6*t_c_hs/x_c_hs)+(100*t_c_hs**4))*((1.34*M**0.18)*(math.cos(Lambda_m)**0.28))
FF_vs = (1+(0.6*t_c_vs/x_c_vs)+(100*t_c_vs**4))*((1.34*M**0.18)*(math.cos(Lambda_m)**0.28))
FF_fuse = 1+(60/f_fuse**3)+(f_fuse/400)
CD_min = ((Cf_wing*FF_wing*Q*S_wet_wing)+(Cf_hs*FF_hs*Q*S_wet_hs)+(Cf_vs*FF_vs*Q*S_wet_vs)+(Cf_fuse*FF_fuse*Q*S_wet_fuse))/S_ref_wing
CD = CD_min+(K*(CL-CL_min)**2)+(CL**2/(math.pi*e*AR))
CD_cruise = CD_min+(K*(CL_cruise-CL_min)**2)+(CL_cruise**2/(math.pi*e*AR))
CD_max = max(CD)

## Plotting
# 3-D Lift
plt.plot(alp, CL) # plots 3-D lift coefficient
plt.grid()
plt.xlabel("Angle of Attack (deg)")
plt.ylabel("3-D Lift Coefficient")
plt.title("3-D Lift Coefficient v. Angle of Attack")
plt.show()

# 3-D Drag
plt.plot(alp, CD) # plots 3-D lift coefficient
plt.grid()
plt.xlabel("Angle of Attack (deg)")
plt.ylabel("3-D Drag Coefficient")
plt.title("3-D Drag Coefficient v. Angle of Attack")
plt.show()

print("Aerodynamic Parameters:")
aero_tab = pt.PrettyTable(["Parameter", "Unit", "Value"])
aero_tab.add_row(["Lift Coefficient (cruise)", "---", CL_cruise])
aero_tab.add_row(["Drag Coefficient (cruise)", "---", CD_cruise])
aero_tab.add_row(["Lift-to-Drag Ratio (cruise)", "---", CL_cruise/CD_cruise])
aero_tab.add_row(["Lift Coefficient (max)", "---", CL_max])
aero_tab.add_row(["Drag Coefficient (max)", "---", CD_max])
aero_tab.add_row(["Lift-to-Drag Ratio (max)", "---", CL_max/CD_max])
aero_tab.add_row(["Lift Coefficient (min)", "---", CL_min])
aero_tab.add_row(["Drag Coefficient (min)", "---", CD_min])
aero_tab.add_row(["Lift-to-Drag Ratio (min)", "---", CL_min/CD_min])

print(aero_tab)

