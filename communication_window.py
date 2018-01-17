#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:42:59 2017

@author: bonar
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import poliastro

from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.plotting import plot
from poliastro.plotting import OrbitPlotter
from poliastro.bodies import Earth, Mars, Sun, Moon

plt.style.use("seaborn")


"""
function which returns viewing angle at a given time
step = time to calculate angle at
T orbital period
r_a apoares
i = inclination
phi_0 = starting position of landing site (0 = same longitude as pericentre)
"""
def visible_angle(step, T, r_a, i, phi_0):
    #astro constants
    T_mars = (24*60+40)*60           #seconds
    G = 6.67408e-11                  #m^3 kg^-1 s^-2
    M_mars = 6.39e23                 #kg
    R_mars = 3390                    #km
    station_lat = 55                 #deg
    mu = 4.282837e13                 #m^3 s^-2
    k = mu*1e-9                      #km^3 s^-2
    
    #orbital elements (orbit defined by apocentre and period)
    r_asc = 37                                        #deg
    arg_p = 90                                        #deg
    true_anom = 0                                     #deg at time t = 0
    a = ((((T**2)/(4*math.pi**2))  * mu  )**(1/3)) * 1e-3  #km
    e = r_a/a-1
    r_p =a*(1-e)
    
    #orbit_classical = [mu, a, e, i, r_asc, arg_p, true_anom]
    #print("a = " +str(a)+ "\ne = "+str(e)+"\ni = " +str(i)+ "\nr_asc = " +str(r_asc)+"\n arg_p = " +str(arg_p)+"\ntrue_nom = " +str(true_anom))
    #print("r_a = "+str(r_a)+"\n r_p = "+str(r_p))
    #T_calc = (2*math.pi*math.sqrt(((a*10**3)**3)/mu))
    
    """
    get position and velocity vector of satellite at t = 0
    then finding vectors at t = step using poliastro twobody propogation
    """
    v_p = math.sqrt(k*((2/(a*(1-e)))-(1/a)))
    r0 = [(r_p * math.cos(math.radians(i))), 0, (r_p * math.sin(math.radians(i)))] * u.km
    v0 = [0, v_p, 0]  * u.km
    
    
    tof = step
    prop = poliastro.twobody.propagation.kepler(k, r0, v0, tof, rtol=1e-10, numiter=35)
    
    r_o = prop[0]
    v_o = prop[1]

  
    """get position vector of landing site at time t, just geometry"""
    w_mars = (2*math.pi)/T_mars    #rad/s
    lander_start = math.degrees(phi_0)
    r_m = R_mars   #km
    theta = math.radians(station_lat)    #inclination from z-axis
    x = r_m * math.sin(theta)*math.cos(phi_0 + w_mars*step)
    y = r_m * math.sin(theta) * math.sin(phi_0 + w_mars*step)
    z = r_m * math.cos(theta)
    r_l = [x, y, z]
   
    """can plot orbits of satellite and an orbit to indicate landing site position"""
    #propped = Orbit.from_vectors(Mars, r_o, v_o*(1/u.s))
    #print('satellite position at time ' + str(step) + ' is ' + str(r_o))
    #op = OrbitPlotter()
    #op.plot(propped)
    #op.plot(Orbit.circular(Mars, alt = 0.1*R_mars * u.km, arglat = (lander_start + math.degrees(w_mars*step)) * u.deg))  
    #plt.savefig(''+(str(step))+'.png')
    #plt.show()

    #print('landing site position at time ' + str(step) + ' is ' + str(r_l))
    
    """use position vectors to find viewing angle""" 
    R_s = []
    for i in range(len(r_o.value)):
        R_s.append(r_o.value[i])        #satellite position
        
    r_dif = []                          #satellite position - lander site position
    for i in range(len(R_s)):
        r_dif.append(R_s[i] - r_l[i])
    
    r_o_mag = np.linalg.norm(R_s)       #magnitude of vectors
    r_l_mag = np.linalg.norm(r_l)
    r_dif_mag = np.linalg.norm(r_dif)
    dot = 0
    
    for i in range(len(r_dif)):         #dot product of vector from station to satellite, and position vector of station (which is normal for sphere)
        dot += (r_dif[i])*(-r_l[i])
    cosX = dot/(r_dif_mag*r_l_mag)


    X = np.arccos(np.clip(cosX, -1.0, 1.0))
    #print(math.degrees(X))

    return((180-math.degrees(X)), a, e, r_p)
    


"""
main loop
"""
#orbit is defined by period, apoares and inclination - other orbital elements are fixed but can be changed in visible_angle function
#can change how orbit is defined by changing top of visible angle function, using apocentre here as I was lowering this
T_mars = (24*60+40)*60   #seconds  
T =  2*T_mars            #orbital period of orbit seconds
r_a = 35000          #apocentre (km)
data_step = 500          #number of measurements per period
cycles = 1               #number of satellite periods to take measurements
inc = -63.435

"""
phi is starting angle of landing site, relative to pericentre of orbit
each degree of phi corresponds to one degree of nodal precession
""" 
for phi in range(0, 360):
    print(str(phi)+"/360")               #deg, inclination from x axis down to perigee in +x direction -z
    phi_0 = math.radians(phi)            #radians, azimuthal angle of landing site from x axis (perigee) at perigee passage 0 is on +x side of mars on x axis
    
    angles = []
    #get viewing angles for every time interval=data_step (angle from normal at station to station-satellite vector)   
    for i in range(data_step*cycles):
        angles.append(visible_angle((T/data_step)*i, T, r_a, inc, phi_0)[0])
        if i ==1:
            a = visible_angle((T/data_step)*i, T, r_a, inc, phi_0)[1]
            e = visible_angle((T/data_step)*i, T, r_a, inc, phi_0)[2]
            r_p = round(visible_angle((T/data_step)*i, T, r_a, inc, phi_0)[3])
            
    """
    print("semimajor axis is " +str(a)+" km")
    print("eccentricity of orbit is " +str(e))
    print("altitude at perigee is " +str(r_p) +'km')
    print("altitude at apogee is " +str(a*(1+math.fabs(e)))+ " km")
    print('period is ' +str(T/60)+ 'minutes')
    """
    x_s = np.arange(data_step*cycles)
    t_s = []
    
    for i in range(len(x_s)):
        t_s.append(x_s[i]/data_step*T/60)
    
    
    #loop through angles and calculate the longest communication window in the first martian day of data
    th_s = [x / 60 for x in t_s]
    windows = [0]
    in_window = False
    for i in range(len(angles)):
        if t_s[i] < T_mars/60:
            if in_window == False:
                if angles[i] > 5 and angles[i]<85:
                    in_window = True
                    start = th_s[i]
            if in_window == True:
                if t_s[(i+2)] > T_mars/60:
                        end = T_mars/60/60
                        wind_length = end-start
                        windows.append(round(wind_length, 4))
                        in_window = False
                elif angles[i] < 5 or angles[i] > 85:
                    in_window = False
                    end = th_s[i]
                    wind_length = end-start
                    windows.append(round(wind_length, 4))
                
                    
    #find longest communication window in first day (sometimes can have more than one window per day)        
    maxwindow = max(windows)

    #plot and show the graph or save figure
    plt.axhline(y=0, color = 'r', linestyle= '-')
    plt.axhline(y=85, color = 'r', linestyle= '-')
    
    plt.plot(th_s, angles)
    for i in range(int(T/T_mars*cycles)):
        plt.axvline(x=(T_mars/60/60)*(i+1), color = 'k', linestyle = ':', alpha = 0.4)
    plt.ylim(-20, 120)
    
    plt.ylabel('Viewing angle (deg)')
    plt.xlabel('Time (hours)')
    plt.title('Visibility of Landing Site\n Changes over time due to nodal precession')
    plt.annotate("Communication window = "+str(maxwindow)+" hours", xy=(0.05, 0.08), xycoords = 'axes fraction', fontsize=16)
    plt.annotate("$r_a = $" +str(r_a)+"km  $r_p = $"+str(r_p)+"km ", xy = (0.05,0.03), xycoords = 'axes fraction', fontsize = 16)
    plt.show()
    """uncomment here to save graphs in same folder as script"""
    #plt.savefig('communication'+str(phi)+'.png')
    plt.clf()
