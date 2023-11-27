# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 16:54:20 2023

@author: Rafael Lima
"""
#%%
import numpy as np
from math import sin
from math import pi
import matplotlib.pyplot as plt
import sympy

#%%
def thomsen_parameters(C11, C13, C33, C55, C66, rhob):
    vp0 = (C33/rho)**(1/2)
    vs0 = (C55/rho)**(1/2)
    epsilon = (C11 - C33)/(2*C33)
    delta = ((C13+C55)**2 - (C33 - C55)**2)/(2*C33*(C33-C55))
    gama = (C66 - C55)/(2*C55)
    
    return vp0, vs0, epsilon, delta, gama


#%%

def degree_to_rad(degree):
    rad = (degree*pi)/360
    return rad

degree = np.linspace(0,90, 1024)
theta = degree_to_rad(degree)
Vp0=3
Vs0=3.5
delta = 0.1
epsilon = 0.2
Vp = np.empty([len(theta)])

Vp0Vs0 = 1.5
rho = 2300
C33 = 25e9
C44 = 20e9
C55 = 12e9
C13 = 16e9


for i in range((len(theta))): 
    Vp[i] = Vp0*(1 + delta*sin(theta[i])**2 + epsilon-delta*sin(theta[i])**4)
    # Vsv[i] = Vs0*(1 + )
#%%
plt.style.use('ggplot')
plt.figure(dpi=300)
plt.subplot(2,2,1)
plt.plot(degree, Vp )
