# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 10:12:59 2022

@author: Giulia
Here is the python code to my project: the modelling of a tiny greenhouse
The goal of this simulation is to determine whether I putting my basil plant 
in a mini greenhouse is enough for it to survive through winter.

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

from dm4bem import read_epw, sol_rad_tilt_surf


l = 0.5             # m length of the cube
Sg = 0.5*0.5       # m² surface of each glass wall
Sc = 5 * Sg    # m² surface of the whole cubic greenhouse (5 faces)
n_nodes = 8
n_branches = 12


wall = {'Conductivity': 1.4,  # W/(m·K)
        'Density': 2500,        # kg/m³
        'Specific heat': 750,  # J/(kg·K)
        'Width': 0.05,
        'Surface': 0.25},  # m²
wall = pd.DataFrame(wall, index=['Glass'])
wall

#Incidence matrix
A = np.zeros([12,8])
A[0, 0] = 1                 # branch 0: -> node 0
A[1, 0], A[1, 5] = -1, 1    # branch 1: node 0 -> node 5
A[2, 1] = 1
A[3, 1], A[3, 5] = -1, 1
A[4, 2] = 1
A[5, 2], A[5, 5] = -1, 1
A[6, 3] = 1
A[7, 3], A[7, 5] = -1, 1
A[8, 4] = 1
A[9, 4], A[9, 5] = -1, 1
A[10, 5], A[10, 6] = -1,1
A[11, 6], A[11, 7] = -1,1

pd.DataFrame(A)



#conductances
G_cd = wall['Conductivity'] / wall['Width'] * wall['Surface']
pd.DataFrame(G_cd, columns={'Conductance'})

#convection
h = pd.DataFrame([{'in': 8., 'out': 25}], index=['h'])  # W/(m²⋅K)
G_cv = h * wall['Surface'][0]     # glass
print(G_cv)

#The equivalence conductance is given by
Ggo = float(1 / (1 / G_cv['out'] + 1 / (2 * G_cd['Glass'])))
Ggi = float(1 / (1 / G_cv['in'] + 1 / (2 * G_cd['Glass'])))

print(Ggo,Ggi)


# We also need to know the conductance and capacity of soil

lambda_Soil = 1 #W/mK
specC_soil = 1480 #J/kgK
h = 0.3
r1 = 0.001
r2 = 0.1
density_soil = 1400 #kg/m3


#conductance in a cylindrical geometry
G_soil = lambda_Soil*np.pi*2*h/(np.log(r2/r1))
C_soil = specC_soil*density_soil*np.pi*r2**2*h
print(G_soil,C_soil)


# Now we can write the conductances matrix $G$

 


G = np.diag([Ggo,Ggi,Ggo,Ggi,Ggo,Ggi,Ggo,Ggi,Ggo,Ggi,int(G_cv['in']),G_soil]) #12x12 matrix
pd.DataFrame(G)


# # Now the capacity matrix

 


C_glass = wall['Density'] * wall['Specific heat'] * wall['Surface'] * wall['Width']                                                       


 


C = np.diag([int(C_glass), int(C_glass), int(C_glass), int(C_glass), int(C_glass), 0, 0, int(C_soil)])
pd.DataFrame(C)


# # Sources: solar ratiation on tilted surfaces

# For each surface we calculate the total irradiation of the sun. To so so we need the orientation of each one of them

 


S_N = {'slope': 90,'azimuth': 180,'latitude': 45} #North
S_E = {'slope': 90,'azimuth': -90,'latitude': 45} #West
S_W = {'slope': 90,'azimuth': 90,'latitude': 45} #East
S_S = {'slope': 90,'azimuth': 0,'latitude': 45} #South
S_Roof = {'slope': 0,'azimuth': 0,'latitude': 45} #Roof


# I download a weather file from https://www.climate.onebuilding.org/

 

filename = './weather_data/CHE_VD_Nyon.Changins.067050_TMYx.2004-2018.epw'
start_date = '2000-01-01'
end_date = '2000-01-05'


 


[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data


 


weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather[(
    weather.index >= start_date) & (
    weather.index < end_date)]
pd.DataFrame(weather)


# Let's try to read the radiation

 


weather[['dir_n_rad', 'dif_h_rad']].plot()
plt.xlabel("Time")
plt.ylabel("Solar radiation (W/m²)")
plt.legend(['$Φ_{direct}$', '$Φ_{diffuse}$'])
plt.show()


# ## Solar radiation for each surface

 


albedo = 0.2
albedo


 


rad_N = sol_rad_tilt_surf(weather, S_N, albedo)
rad_W = sol_rad_tilt_surf(weather, S_W, albedo)
rad_S = sol_rad_tilt_surf(weather, S_S, albedo)
rad_E = sol_rad_tilt_surf(weather, S_E, albedo)
rad_Roof = sol_rad_tilt_surf(weather, S_Roof, albedo)

rad_S.plot()
plt.xlabel("Time")
plt.ylabel("Solar radiation (W/m²)")
plt.show()


# We have the values in $W/m^2$, but we want the values in $[W]$.
# To do so, we must do $\Phi_s = \alpha * Surface * \Phi$
# 
# Where the surface is the same for each surface (=Sg)

# The total radiation is the sum of direct diffuse and reflected

 


rad_N['EtotN'] = rad_N.sum(axis=1)
rad_S['EtotS'] = rad_S.sum(axis=1)
rad_E['EtotE'] = rad_E.sum(axis=1)
rad_W['EtotW'] = rad_W.sum(axis=1)
rad_Roof['EtotR'] = rad_Roof.sum(axis=1)


 


rad_N.plot()


# 

 





# ## Now that we have the radiation on the 5 surfaces, we can set up the flow and temperature source vectors, b and f
# 
# we set the values at 1 where there is a temperature or a flow source. 

 


b = np.zeros(12)        # branches
b[[0, 2, 4,6,8]] = 1   # branches with temperature sources
print(f'b = ', b)

f = np.zeros(8)         # nodes
f[[0, 1, 2, 3, 4, 6]] = 1     # nodes with heat-flow sources
print(f'f = ', f)

y = np.zeros(8)         # nodes
y[[5]] = 1              # nodes (temperatures) of interest
print(f'y = ', y)


# State space representation

[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)

# Steady-tate analysis, to check that our model is not wrong
# We set the outdoor temperature source to 10°C and the flows to zero. 
# In the steady state we expect all the temperatures to be at the outdoor temperature

b[[0, 2, 4 ,6 ,8]] = 10     # outdoor temperature          
f = np.zeros(8)            # flow-rate sources
θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
print(f'θ = {θ} °C')


# Now we can set up the vector u, with the non-zero elements of b and f.

bT = np.array([10, 10, 10, 10, 10]) #To, To, To, Tisp]
fQ = np.array([0, 0, 0, 0, 0, 0])         
u = np.hstack([bT, fQ])
print(f'u = {u}')
yss = (-Cs @ np.linalg.inv(As) @ Bs + Ds) @ u
print(f'yss = {yss} °C')


# Dynamic simulation


λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
print('Time constants: \n', -1 / λ, 's \n')
print('2 x Time constants: \n', -2 / λ, 's \n')
dtmax = min(-2. / λ)
print(f'Maximum time step: {dtmax:.2f} s = {dtmax / 60:.2f} min')


dt=round((dtmax)/1.2868)
print(f'The time step we use for our simulation is: {dt} s')

#settling time
t_resp = 4 * max(-1 / λ)
print('Time constants: \n', -1 / λ, 's \n')
print(f'Settling time: {t_resp:.0f} s = {t_resp / 60:.1f} min = {t_resp / (3600):.2f} h = {t_resp / (3600 * 24):.2f} days')

# Step response

duration = round(t_resp)            # seconds, larger than response time
n = int(np.floor(duration / dt))    # number of time steps
t = np.arange(0, n * dt, dt)        # time vector for n time steps

print(f'Duration = {duration} s')
print(f'Number of time steps = {n}')
pd.DataFrame(t, columns=['time'])


#step response in the steady state analysis To=O and flows are zeros

u = np.zeros([11, n])                # u = [To To To Tisp Φo Φi Qa Φa
u[0:5, :] = 10 * np.ones([5, n])    # To = 10 for n time steps
print('u = ')
pd.DataFrame(u)

n_s = As.shape[0]                      # number of state variables
θ_exp = np.zeros([n_s, t.shape[0]])    # explicit Euler in time t
θ_imp = np.zeros([n_s, t.shape[0]])    # implicit Euler in time t

I = np.eye(n_s)                        # identity matrix

for k in range(n - 1):
    θ_exp[:, k + 1] = (I + dt * As) @        θ_exp[:, k] + dt * Bs @ u[:, k]
    θ_imp[:, k + 1] = np.linalg.inv(I - dt * As) @        (θ_imp[:, k] + dt * Bs @ u[:, k])


y_exp = Cs @ θ_exp + Ds @  u
y_imp = Cs @ θ_imp + Ds @  u


fig, ax = plt.subplots()
ax.plot(t / 3600, y_exp.T, t / 3600, y_imp.T)
ax.set(xlabel='Time [h]',
       ylabel='$T_i$ [°C]',
       title='Step input: To')
ax.legend(['Implicit', 'Explicit'])
plt.show()


 


data = pd.concat([weather['temp_air'], rad_N['EtotN'], rad_E['EtotE'], rad_S['EtotS'],rad_W['EtotW'],rad_Roof['EtotR']], axis=1)
weather['temp_air'].plot()

#resampling
data = data.resample(str(dt) + 'S').interpolate()
data['temp_air'].plot()


pd.DataFrame(data)


To = data['temp_air']
#absorptivity of glass

α_g = 0.38    # short wave absortivity: reflective blue glass
τ_g = 0.30    # short wave transmitance: reflective blue glass

Φ_N = data['EtotN']*Sg*α_g
Φ_W = data['EtotW']*Sg*α_g
Φ_E = data['EtotE']*Sg*α_g
Φ_S = data['EtotS']*Sg*α_g
Φ_Roof = data['EtotR']*Sg*α_g
Φ_pot = (data['EtotN']+data['EtotE']+data['EtotS']+data['EtotW']+data['EtotR'])*Sg*τ_g

u = pd.concat([To, To, To, To, To, Φ_N, Φ_W, Φ_S, Φ_E, Φ_Roof, Φ_pot], axis=1)
u.columns.values[[0,1,2,3,4]] = ['To','To','To','To','To']
u.columns.values[[5,6,7,8,9,10]] = ['$\Phi_N$','$\Phi_W$','$\Phi_S$','$\Phi_E$','$\Phi_R$','$\Phi_pot$']

pd.DataFrame(u)


# Initial conditions

θ_exp = 10 * np.ones([As.shape[0], u.shape[0]])

for k in range(u.shape[0] - 1):
    θ_exp[:, k + 1] = (I + dt * As) @ θ_exp[:, k]        + dt * Bs @ u.iloc[k, :]



y_exp = Cs @ θ_exp + Ds @ u.to_numpy().T

t = dt * np.arange(data.shape[0])   # time vector

fig, axs = plt.subplots(2, 1)
# plot indoor and outdoor temperature
axs[0].plot(t / 3600 / 24, y_exp[0, :], label='$T_{indoor}$')
axs[0].plot(t / 3600 / 24, data['temp_air'], label='$T_{outdoor}$')
axs[0].set(xlabel='Time [days]',
           ylabel='Temperatures [°C]',
           title='Simulation for weather')
axs[0].legend(loc='upper right')

# plot total solar radiation 
axs[1].plot(t / 3600 / 24, data['EtotN'], label='$Φ_{total}$')
axs[1].set(xlabel='Time [days]',
           ylabel='Heat flows [W]')
axs[1].legend(loc='upper right')

fig.tight_layout()