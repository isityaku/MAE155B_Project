### MAE155B Sizing and Scoring Analysis
# Made By: Lance De La Cruz

## Outline
# 1. Read .csv file of weight distribution
# 2. Calculate Center of Gravity, Neutral Point, Static Margin
# 3. Calculate Aerodynamic Center (AC)
# 4. Calculate Tail Lift-Curve Slope and Coefficients

## Initializations
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import prettytable as pt

os.system('cls') # Clearing terminal

## CSV Reading
w_data = pd.read_csv("MAE155B_Stability/MAE155B_Weight_Breakdown.csv")
w_data_col = w_data.columns.to_numpy()
w_data_name = w_data["Component Name"].to_numpy()
w_data_w = w_data["Weight [g]"].to_numpy()
w_data_x = w_data["X [m]"].to_numpy()
w_data_y = w_data["Y [m]"].to_numpy()
w_data_z = w_data["Z [m]"].to_numpy()
w_data_mat = w_data["Material"].to_numpy()

## Moment Calculation
wx = []
wy = []
wz = []
for i in range(len(w_data_name)):
    wx.append(w_data_w[i]*w_data_x[i])
    wy.append(w_data_w[i]*w_data_y[i])
    wz.append(w_data_w[i]*w_data_z[i])
Ewi = np.sum(w_data_w)
Ewxi = np.sum(wx)
Ewyi = np.sum(wy)
Ewzi = np.sum(wz)

print(Ewxi)
