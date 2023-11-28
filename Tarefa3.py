# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 12:25:09 2023

@author: Rafael Lima
"""

def gassman(Vp_dry, Vs_dry, rho_dry, rho_fl_sat, K_fl_sat, K_min, phi):
    import numpy as np
    rho_sat = rho_dry + phi*rho_fl_sat
    
    mi_dry = rho_dry*Vs_dry**2
    K_dry = rho_dry*Vp_dry**2-(4/3)*mi_dry
    K_dry = K_dry/1e9
    K_sat = K_dry +( ((1-(K_dry/K_min))**2)/(phi/K_fl_sat + (1-phi)/K_min + K_dry/K_min**2))
    
    mi_sat = mi_dry
    
    Vp_sat = np.sqrt(((K_sat + (4/3)*mi_sat)/rho_sat))
    Vs_sat = np.sqrt((mi_sat/rho_sat))
    

    return Vp_sat, Vs_sat, rho_sat, K_sat, K_dry
#%% 

import pandas as pd

from io import StringIO
import numpy as np
# Define the table as a string
table_data = """
Sample,Density Dry,Vp(Radial) Dry,Vs(Radial) Dry,Porosity,Permeabilidade (mD),Density Sat,Vp (Radial) Sat,Vs (Radial) Sat
AB2,2214.162835,3721.618954,2365.119197,20.15,5.96,2422.546207,3850.868233,2376.8923
AB4,2128.969011,3459.633028,2291.008505,23.42,6.17,2364.802857,3657.613967,2303.78659
BC1,2239.489443,4068.965517,2463.144162,18.94,4.82,2424.505795,4154.015402,2495.3478
BC2,2209.519454,3831.129196,2404.853129,19.42,5.72,2414.237486,3908.260689,2420.6545
BC4,2120.34854,3468.26127,2472.131148,22.92,14.9,2358.81981,3584.672435,2600
BC8,2125.229544,3477.460902,2139.219015,23.07,5.64,2356.833365,3589.74359,2158.2
BC12,2271.428403,4370.37037,2539.340955,18.34,1.18,2447.95582,4416.374269,2554.803789
CD9,2064.000631,3314.360771,2070.568928,24.24,4.47,2317.481052,3428.442029,2104.45987
EF02,2002.204663,3282.722513,2109.927089,25.08,1.74,2255.019028,3476.89464,2200
"""

# Use StringIO to simulate reading from a file-like object
table_file = StringIO(table_data)

# Read the table into a Pandas DataFrame
df = pd.read_csv(table_file)

# Display the DataFrame
print(df)


#%%



# rho_dry = 2214.162835
rho_dry = np.array([df['Density Dry']]).reshape(-1, 1)
rho_fl_sat = 997 # Densidade da água em kg/m³
K_fl_sat = 2.05 #GPa modulo de incompressibilidade da água
phi = np.array(df['Porosity']/100).reshape(-1,1)
rho_rock = rho_dry/(1-phi)
# Vp_dry = np.array(df['Vp(Radial) Dry']).reshape(-1,1)
# Vs_dry = np.array(df['Vs(Radial) Dry']).reshape(-1,1)

# # Vp_sat = np.array(df[])

Vp_sat_exp = np.array(df['Vp (Radial) Sat']).reshape(-1,1)
Vp_dry_exp = np.array(df['Vp(Radial) Dry']).reshape(-1,1)
Vs_sat_exp = np.array(df['Vs (Radial) Sat']).reshape(-1,1)
Vs_dry_exp = np.array(df['Vs(Radial) Dry']).reshape(-1,1)

rho_sat_exp = np.array(df['Density Sat']).reshape(-1,1)
K_sat_exp = (Vp_sat_exp**2 - (4/3)*Vs_sat_exp**2)*rho_sat_exp
K_sat_exp = K_sat_exp/1e9
#%%
# phi = 0.215
# phi = phi/100

rho_1 = 3100 # Ankerita
rho_2 = 2710 # Calcita
# rho_3 = 2870
# rho_medio = (rho_1+ rho_2+rho_3)/3
# rho_rock = 2214.162835

K_1 = 70
K_2 = 76
# K_3 = 94

mi_1 = 32
mi_2 = 30
mi_3 = 45

f_1 = (rho_rock - rho_2)/(rho_1 - rho_2)
f_2 = 1 - f_1


# Vp_dry = 3721.618954
# Vs_dry = 2365.119197

K_min = f_1*K_1 + f_2*K_2 
mi_min = f_1*mi_1 + f_2*mi_2

# K_eff = (Vp**2 - 4/3*Vs**2)*rho_rock

#%%

Vp_sat_gass, Vs_sat_gass, rho_sat_gass, K_sat_gass, K_dry = gassman(Vp_dry_exp, Vs_dry_exp, rho_dry, rho_fl_sat, K_fl_sat, K_min, phi)

n = np.arange(1,len(Vp_dry_exp)+1)
#%%

import matplotlib.pyplot as plt

plt.figure(dpi=300, figsize=[7,7])
plt.style.use('ggplot')
plt.subplot(221)
plt.scatter(n, Vp_dry_exp, label = 'Dry Experimental')
plt.scatter(n, Vp_sat_exp, label = 'Sat Experimental')
plt.scatter(n, Vp_sat_gass, label = 'Sat Gassman')
plt.ylabel('Velocity[$m/s$]')
plt.xlabel('Samples')
plt.legend(fontsize = 7)
plt.title('$V_p$')
plt.subplot(222)
plt.scatter(n, Vs_dry_exp, label = 'Dry Experimental')
plt.scatter(n, Vs_sat_exp, label = 'Sat Experimental')
plt.scatter(n, Vs_sat_gass, label = 'Sat Gassman')
plt.ylabel('Velocity[$m/s$]')
plt.xlabel('Samples')
plt.legend(fontsize = 7, loc='best')
plt.title('$V_s$')
plt.subplot(223)
plt.scatter(n, rho_dry, label = 'Dry Experimental')
plt.scatter(n, rho_sat_exp, label = 'Sat Experimental')
plt.scatter(n, rho_sat_gass, label = 'Sat Gassman')
plt.ylabel('Density[$kg/m^3]$')
plt.xlabel('Samples')
plt.legend(fontsize = 7)
plt.title('$\u03c1$')
plt.subplot(224)
plt.scatter(n, K_dry, label = 'Dry Experimental')
plt.scatter(n, K_sat_exp, label = 'Sat Experimental')
plt.scatter(n, K_sat_gass, label = 'SatGassman')
plt.ylabel('Pressure[$GPa$]')
plt.xlabel('Samples')
plt.legend(fontsize = 7)
plt.title('$K_{sat}$')
plt.suptitle('Water Saturation - Gassman', fontsize=20)
plt.tight_layout()


#%% Voitgh e Reuss


plt.figure(dpi=300, figsize=[8,8])
plt.style.use('ggplot')
plt.subplot(221)
plt.scatter(phi, Vp_dry_exp, label = 'Dry Experimental')
plt.scatter(phi, Vp_sat_exp, label = 'Sat Experimental')
plt.scatter(phi, Vp_sat_gass, label = 'Sat Gassman')
plt.ylabel('Velocity[$m/s$]')
plt.xlabel('Porosity[$v/v$]')
plt.legend(fontsize = 7)
plt.title('$V_p$')
plt.subplot(222)
plt.scatter(phi, Vs_dry_exp, label = 'Dry Experimental')
plt.scatter(phi, Vs_sat_exp, label = 'Sat Experimental')
plt.scatter(phi, Vs_sat_gass, label = 'Sat Gassman')
plt.ylabel('Velocity[$m/s$]')
plt.xlabel('Porosity[$v/v$]')
plt.legend(fontsize = 7, loc='best')
plt.title('$V_s$')
plt.subplot(223)
plt.scatter(phi, rho_dry, label = 'Dry Experimental')
plt.scatter(phi, rho_sat_exp, label = 'Sat Experimental')
plt.scatter(phi, rho_sat_gass, label = 'Sat Gassman')
plt.ylabel('Density[$kg/m^3]$')
plt.xlabel('Porosity[$v/v$]')
plt.legend(fontsize = 7)
plt.title('$\u03c1$')
plt.subplot(224)
plt.scatter(phi, K_dry, label = 'Dry Experimental')
plt.scatter(phi, K_sat_exp, label = 'Sat Experimental')
plt.scatter(phi, K_sat_gass, label = 'SatGassman')
plt.ylabel('Pressure[$GPa$]')
plt.xlabel('Porosity[$v/v$]')
plt.legend(fontsize = 7)
plt.title('$K_{sat}$')
plt.suptitle('Water Saturation - Gassman', fontsize=20)
plt.tight_layout()
#%% Voitgh e Reuss

space_porosity = np.arange(0,1.1,0.1)
K_reuss = np.array([])
# for i in range(0, len(K_min)):
#     K_reuss = np.append(K_reuss, phi[i]/K_min[i])
#     K_reuss = 1/K_reuss
# K_reuss = 1/K_reuss

K_reuss = ((1 - space_porosity)/K_min[0] + (space_porosity)/K_fl_sat)**-1
K_voigt = ((1 - space_porosity)*K_min[0] + space_porosity*K_fl_sat)

# K_voigt = (K_1 * K_2) * (1 - space_porosity) / (K_1 - K_2) 

plt.figure(dpi=300)
plt.plot(space_porosity, K_reuss, label='Reuss')
plt.plot(space_porosity, K_voigt, label='Voigt')
plt.scatter(phi, K_dry, label = 'Dry Experimental')
plt.scatter(phi, K_sat_exp, label = 'Sat Experimental')
plt.scatter(phi, K_sat_gass)

