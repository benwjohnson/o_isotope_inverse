#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:13:28 2019

@author: benjohnson
"""

#### Revised inversion compilation for Johnson and Wing, 2019. Herein, we invert 
#O-isotope data from a number of altered ocean crust sections as well as forward-
# modelled results to demonstrate the validity of our system-wide mass balance approach



import pandas as pd 
import numpy as np
import scipy as sp

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset, zoomed_inset_axes)

from o_isotope_invert import o_isotope_invert

plt.close('all')

font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


d18Oinit_morb = 5.5
#%% Calibrating on  forward models 
d18Oinit_hoku = 7.5


#Kuroko (Cathles, 1983) 
hoku_path = 'hokuroku_50kyr.csv'
data_raw_hoku50 = pd.read_csv(hoku_path)
del18O_hoku50 = pd.Series(data_raw_hoku50['d18O']).values
del18O_hoku50 = np.add(del18O_hoku50,2)
x_hoku50 = pd.Series(data_raw_hoku50['x']).values
z_hoku50 = pd.Series(data_raw_hoku50['z']).values
temp_hoku50 = pd.Series(data_raw_hoku50['temp']).values
temp_hoku50 = np.multiply(temp_hoku50,0.85)


d18Oinit_hoku50 = d18Oinit_hoku
num_runshoku50 =  del18O_hoku50.size

  

hoku50_water_initial = np.zeros(num_runshoku50)
hoku50_water_outgoing = np.zeros(num_runshoku50)
hoku50_W_R = np.zeros(num_runshoku50)
hoku50_final_rock = np.zeros(num_runshoku50)
hoku50_moles_fluid = np.zeros(num_runshoku50)
hoku50_moles_Orock = np.zeros(num_runshoku50)
for irun in range(0,num_runshoku50):
    del18O_in_hoku50 = del18O_hoku50
    Temp_in_hoku50 = temp_hoku50
    x_hoku50_in = x_hoku50
    z_hoku50_in = z_hoku50
    
    num_takeout = np.random.randint(0,len(del18O_in_hoku50))
    
    del18O_in_hoku50 = np.delete(del18O_in_hoku50,num_takeout)
    Temp_in_hoku50 =np.delete(Temp_in_hoku50,num_takeout)
    x_hoku50_in =  np.delete(x_hoku50,num_takeout)
    z_hoku50_in =np.delete(z_hoku50,num_takeout)
    
    [hoku50_iso_grid,hoku50_x_grid,hoku50_z_grid,conc_matrix_hoku50,change_matrix_hoku50,qsolved_hoku50,\
     moles_fluid_hoku50,W_R_taper,water_initial_hoku50, water_outgoing_hoku50, temp_test_hoku50, temp_grid_hoku50,moles_Orock_hoku50,hoku50_epsilon_r_w] = \
    o_isotope_invert(del18O_in_hoku50,Temp_in_hoku50,x_hoku50_in,z_hoku50_in,d18Oinit_hoku50)
    hoku50_water_initial[irun] = water_initial_hoku50
    hoku50_water_outgoing[irun] = water_outgoing_hoku50
    hoku50_moles_fluid[irun] = moles_fluid_hoku50
    hoku50_W_R[irun] = W_R_taper
    hoku50_final_rock[irun] = del18O_in_hoku50.mean()
    hoku50_moles_Orock[irun] = moles_Orock_hoku50

# Wing and Ferry (2007)
wing_path = 'wing_and_ferry.csv'   
data_raw_wing = pd.read_csv(wing_path)
del18O_wing = pd.Series(data_raw_wing['d18O']).values
#del18O_wing = np.multiply(del18O_wing,1.01)
x_wing = pd.Series(data_raw_wing['x']).values
#x_wing = np.multiply(x_wing, 1000)
z_wing = pd.Series(data_raw_wing['z']).values
temp_wing = pd.Series(data_raw_wing['temp']).values

num_runswing = del18O_wing.size



d18Oinit_wing = 5.7 #in Wing and Ferry, it's listed as 5.5permil, but Boz found his code and it's actually 5.7



wing_water_initial = np.zeros(num_runswing)
wing_water_outgoing = np.zeros(num_runswing)
wing_moles_fluid = np.zeros(num_runswing)
wing_W_R = np.zeros(num_runswing)
wing_final_rock = np.zeros(num_runswing)
wing_moles_Orock = np.zeros(num_runswing)
for irun in range(0,num_runswing):
    del18O_in_wing = del18O_wing
    Temp_in_wing = temp_wing
    x_wing_in = x_wing
    z_wing_in = z_wing
    
    num_takeout = np.random.randint(0,len(del18O_in_wing))
   
    del18O_in_wing = np.delete(del18O_in_wing,num_takeout)
    Temp_in_wing =np.delete(Temp_in_wing,num_takeout)
    x_wing_in =  np.delete(x_wing,num_takeout)
    z_wing_in =np.delete(z_wing,num_takeout)
    
    [wing_iso_grid,wing_x_grid,wing_z_grid,conc_matrix_wing,change_matrix_wing,qsolved_wing,
     moles_fluid_wing,W_R_taper,water_initial_wing, water_outgoing_wing, temp_test_wing, temp_grid_wing,moles_Orock_wing,wing_epsilon_r_w] = \
    o_isotope_invert(del18O_in_wing,Temp_in_wing,x_wing_in,z_wing_in,d18Oinit_wing)
    wing_water_initial[irun] = water_initial_wing
    wing_water_outgoing[irun] = water_outgoing_wing
    wing_moles_fluid[irun] = moles_fluid_wing
    wing_W_R[irun] = W_R_taper
    wing_final_rock[irun] = del18O_in_wing.mean()
    wing_moles_Orock[irun] = moles_Orock_wing


#Skaergaard
skaer_path = 'skaergaard_model.csv'
data_raw_skaer = pd.read_csv(skaer_path)
del18O_skaer = pd.Series(data_raw_skaer['d18O']).values
x_skaer = pd.Series(data_raw_skaer['x']).values
x_skaer = np.multiply(x_skaer, 1000)
z_skaer = pd.Series(data_raw_skaer['z']).values
temp_skaer = pd.Series(data_raw_skaer['temp']).values

num_runsskaer = del18O_skaer.size

temp_skaer = np.multiply(temp_skaer,0.91)
#

d18Oinit_skaer = 7
  
skaer_water_initial = np.zeros(num_runsskaer)
skaer_water_outgoing = np.zeros(num_runsskaer)
skaer_moles_fluid = np.zeros(num_runsskaer)
skaer_W_R = np.zeros(num_runsskaer)
skaer_final_rock = np.zeros(num_runsskaer)
skaer_moles_Orock = np.zeros(num_runsskaer)

for irun in range(0,num_runsskaer):
    del18O_in_skaer = del18O_skaer
    Temp_in_skaer = temp_skaer
    x_skaer_in = x_skaer
    z_skaer_in = z_skaer
    
    num_takeout = np.random.randint(0,len(del18O_in_skaer))
    
    del18O_in_skaer = np.delete(del18O_in_skaer,num_takeout)
    Temp_in_skaer =np.delete(Temp_in_skaer,num_takeout)
    x_skaer_in =  np.delete(x_skaer,num_takeout)
    z_skaer_in =np.delete(z_skaer,num_takeout)
    
    [skaer_iso_grid,skaer_x_grid,skaer_z_grid,conc_matrix_skaer,change_matrix_skaer,qsolved_skaer,
     moles_fluid_skaer,W_R_taper,water_initial_skaer, water_outgoing_skaer, temp_test_skaer, temp_grid_skaer,moles_Orock_skaer,skaer_epsilon_r_w] = \
    o_isotope_invert(del18O_in_skaer,Temp_in_skaer,x_skaer_in,z_skaer_in,d18Oinit_skaer)
    skaer_water_initial[irun] = water_initial_skaer
    skaer_water_outgoing[irun] = water_outgoing_skaer
    skaer_moles_fluid[irun] = moles_fluid_skaer
    skaer_W_R[irun] = W_R_taper
    skaer_final_rock[irun] = del18O_in_skaer.mean()    
    skaer_moles_Orock[irun] = moles_Orock_skaer
#%% Inverting modern/young sections
    
#East Pacific Rise
data_raw = pd.read_csv('/Users/benjohnson/Desktop/CU_google_drive/current_projects/geochemical_inverse_modeling/EPR/gillis_et_al_data.csv')
del18Oepr = pd.Series(data_raw['d18O']).values

x = pd.Series(data_raw['x']).values
x = np.divide(x,1e2)
y = pd.Series(data_raw['y']).values
y = np.divide(y,1e2)
temp = np.multiply(np.divide(1,del18Oepr**(0.5)),700)   

num_runsEPR = del18Oepr.size

#get height above lowest sample from x/y data and slope
slope = 52 #slope of outcrop in Gillis et al., 2001 Fig 2 in degrees
intercept = 1900 #m away from SW corne of plot area in Gillis et al., 2001 Fig2
thickness = 1710 #total thickness
distance = (np.subtract(intercept,y)**2)**0.5 #distance from each sample point to plane of projection
dft = np.multiply(distance,np.tan(slope)) #distance from the top of the section
height = np.subtract(thickness,dft) 



EPR_water_initial = np.zeros(num_runsEPR)
EPR_water_outgoing = np.zeros(num_runsEPR)
EPR_moles_fluid = np.zeros(num_runsEPR)
EPR_W_R = np.zeros(num_runsEPR)
EPR_final_rock = np.zeros(num_runsEPR)
EPR_moles_Orock = np.zeros(num_runsEPR)
for irun in range(0,num_runsEPR):
    del18O_in_EPR = del18Oepr
    Temp_in_EPR = temp
    x_EPR = x
    z_EPR = height
    
    num_takeout = np.random.randint(0,len(del18O_in_EPR))

#    
    del18O_in_EPR = np.delete(del18O_in_EPR,num_takeout)
    Temp_in_EPR =np.delete(Temp_in_EPR,num_takeout)
    x_EPR =  np.delete(x_EPR,num_takeout)
    z_EPR =np.delete(z_EPR,num_takeout)
    
    [EPR_iso_grid,EPR_x_grid,EPR_z_grid,conc_matrix_EPR,change_matrix_EPR,
     qsolved_EPR,moles_fluid_EPR,W_R_taper,water_initial_EPR, water_outgoing_EPR,temp_test_EPR,temp_grid_EPR,moles_Orock_EPR, EPR_epsilon_r_w] = \
    o_isotope_invert(del18O_in_EPR,Temp_in_EPR,x_EPR,z_EPR,d18Oinit_morb)
    EPR_water_initial[irun] = water_initial_EPR
    EPR_water_outgoing[irun] = water_outgoing_EPR
    EPR_moles_fluid[irun] = moles_fluid_EPR
    EPR_W_R[irun] = W_R_taper
    EPR_final_rock[irun] = del18O_in_EPR.mean()
    EPR_moles_Orock[irun] = moles_Orock_EPR
EPR_water_moles = np.linalg.norm(qsolved_EPR)
#Hokoroku
ohmoto_path = 'green_et_al.csv'
data_raw_ohmoto = pd.read_csv(ohmoto_path)
del18O_ohmoto = pd.Series(data_raw_ohmoto['d18O']).values
#del18O_ohmoto = np.add(del18O_ohmoto,-1)
x_ohmoto = pd.Series(data_raw_ohmoto['x']).values
z_ohmoto = pd.Series(data_raw_ohmoto['z']).values
temp_ohmoto = pd.Series(data_raw_ohmoto['temp']).values
temp_ohmoto = np.multiply(temp_ohmoto,0.85)
#
#d18Oinit_hoku = d18Oinit_morb   
d18Oinit_ohmoto = d18Oinit_hoku
num_runsohmoto = del18O_ohmoto.size
  
ohmoto_water_initial = np.zeros(num_runsohmoto)
ohmoto_water_outgoing = np.zeros(num_runsohmoto)
ohmoto_moles_fluid = np.zeros(num_runsohmoto)
ohmoto_W_R = np.zeros(num_runsohmoto)
ohmoto_final_rock = np.zeros(num_runsohmoto)
ohmoto_moles_Orock = np.zeros(num_runsohmoto)

for irun in range(0,num_runsohmoto):
    del18O_in_ohmoto = del18O_ohmoto
    Temp_in_ohmoto = temp_ohmoto
    x_ohmoto_in = x_ohmoto
    z_ohmoto_in = z_ohmoto
    
    num_takeout = np.random.randint(0,len(del18O_in_ohmoto))
    
    del18O_in_ohmoto = np.delete(del18O_in_ohmoto,num_takeout)
    Temp_in_ohmoto =np.delete(Temp_in_ohmoto,num_takeout)
    x_ohmoto_in =  np.delete(x_ohmoto,num_takeout)
    z_ohmoto_in =np.delete(z_ohmoto,num_takeout)
    
    [ohmoto_iso_grid,ohmoto_x_grid,ohmoto_z_grid,conc_matrix_ohmoto,change_matrix_ohmoto,
     qsolved_ohmoto,moles_fluid_ohmoto,W_R_taper,water_initial_ohmoto, water_outgoing_ohmoto, temp_test_ohmoto, temp_grid_ohmoto,moles_Orock_ohmoto, ohmoto_epsilon_r_w] = \
    o_isotope_invert(del18O_in_ohmoto,Temp_in_ohmoto,x_ohmoto_in,z_ohmoto_in,d18Oinit_ohmoto)
    ohmoto_water_initial[irun] = water_initial_ohmoto
    ohmoto_water_outgoing[irun] = water_outgoing_ohmoto
    ohmoto_moles_fluid[irun] = moles_fluid_ohmoto
    ohmoto_W_R[irun] = W_R_taper
    ohmoto_final_rock[irun] = del18O_in_ohmoto.mean()
    ohmoto_moles_Orock[irun] = moles_Orock_ohmoto
  

#Solea Graben
solea_path = 'solea_graben.csv'
data_raw_solea = pd.read_csv(solea_path)
del18O_solea = pd.Series(data_raw_solea['d18O']).values
x_solea = pd.Series(data_raw_solea['x']).values
x_solea = np.multiply(x_solea, 1000)
z_solea = pd.Series(data_raw_solea['z']).values
temp_solea = pd.Series(data_raw_solea['temp']).values


temp_solea = np.multiply(temp_solea,0.88)
#

d18Oinit_solea = d18Oinit_morb
num_runssolea = del18O_solea.size  

solea_water_initial = np.zeros(num_runssolea)
solea_water_outgoing = np.zeros(num_runssolea)
solea_moles_fluid = np.zeros(num_runssolea)
solea_W_R = np.zeros(num_runssolea)
solea_final_rock = np.zeros(num_runssolea)
solea_moles_Orock = np.zeros(num_runssolea)

for irun in range(0,num_runssolea):
    del18O_in_solea = del18O_solea
    Temp_in_solea = temp_solea
    x_solea_in = x_solea
    z_solea_in = z_solea
    
    num_takeout = np.random.randint(0,len(del18O_in_solea))
  
    del18O_in_solea = np.delete(del18O_in_solea,num_takeout)
    Temp_in_solea =np.delete(Temp_in_solea,num_takeout)
    x_solea_in =  np.delete(x_solea,num_takeout)
    z_solea_in =np.delete(z_solea,num_takeout)
    
    [solea_iso_grid,solea_x_grid,solea_z_grid,conc_matrix_solea,change_matrix_solea,
     qsolved_solea,moles_fluid_solea,W_R_taper,water_initial_solea, water_outgoing_solea,
     temp_test_solea, temp_grid_solea,moles_Orock_solea,solea_epsilon_r_w] = \
    o_isotope_invert(del18O_in_solea,Temp_in_solea,x_solea_in,z_solea_in,d18Oinit_solea)
    solea_water_initial[irun] = water_initial_solea
    solea_water_outgoing[irun] = water_outgoing_solea
    solea_moles_fluid[irun] = moles_fluid_solea
    solea_W_R[irun] = W_R_taper
    solea_final_rock[irun] = del18O_in_solea.mean()
    solea_moles_Orock[irun] = moles_Orock_solea
solea_water_moles = np.linalg.norm(qsolved_solea)

#%% Panorama
geochem = pd.read_csv('/Users/benjohnson/Desktop/CU_google_drive/current_projects/geochemical_inverse_modeling/panorama_vms_geochemistry/practice_python/panorama_majors_oxygen.csv')


#spatial coordinates
x_pan=pd.Series(geochem['mE']).values
y_pan=pd.Series(geochem['mN']).values; 

#O-isotopes and temperature
del18Opan = pd.Series(geochem['d18O']).values
temppan = pd.Series(geochem['Temp']).values

#initial d18O value
d18O_pano = d18Oinit_morb
num_runspano = del18Opan.size
##define extent of each cell, then patch together to make whole panorama O isotope and temp contour 

#Kangaroo Caves
x_redux = np.zeros(len(x_pan))
y_redux = np.zeros(len(y_pan)) 

for i in range(len(x_pan)):
    if x_pan[i]>729000 and x_pan[i]<733000 and y_pan[i]>7650000 and y_pan[i]<7655000:
        x_redux[i] = x_pan[i]
        y_redux[i] = y_pan[i]
    #remove zeros

indices=np.nonzero(x_redux)[0]

x_redux=x_redux[x_redux!=0]
y_redux=y_redux[y_redux!=0]  

O_isotope_redux=del18Opan[indices]#geochem.d18O[indices] 
Temp_redux=temppan[indices] #geochem.Temp[indices]
#Correct for dip, project points onto plane, essentially making a true cross section to grid
#First, define a line for the base of the section (md), just based off of measuring on the map
dip_angle=55 #dip of beds

A=[7654466,730016]; B=[7649666,730801] #First entry is northing, second easting
md=(A[0]-B[0])/(A[1]-B[1])#slope of the line
bd=A[0]-md*A[1] #y-intercept

#then, for every point, find its position (D) along that line, by defining a line P that is perpendicular
# to D
mp=-1/md #slope of this line is the same for all points
bp=np.subtract(y_redux,mp*x_redux) #y-intercept is unique for each point
#then, at each projection of point onto base line D, line D = line P
PosDx = (bd-bp)/(mp-md) #x-position in original axes
PosDy = mp*PosDx + bp #y-position in original axes
PosDtrue = ((A[0]-PosDx)**2+(A[1]-PosDy)**2)**0.5#position along the line at the base of the formation, using A as the origin in the NW

#then, need to find the vertical distance above the line D for data points projected perpendicularly
tempx = np.subtract(PosDx,x_redux)**2
tempy = np.subtract(PosDy,y_redux)**2
DistOrig=np.add(tempx,tempy)**0.5 #Distance from point on line D to original point
 #calculate distance from line D to point Z, which is stratigraphic height
Z=np.tan(np.deg2rad(dip_angle))*DistOrig
Z=np.subtract(Z,Z.min())
#--Sulphur Springs

x_redux_SS = np.zeros(len(x_pan))
y_redux_SS = np.zeros(len(y_pan)) 
for i in range(len(x_pan)):
    if x_pan[i]>727000 and x_pan[i]<735000 and y_pan[i]>7657000 and y_pan[i]<7660000:
        x_redux_SS[i] = x_pan[i]
        y_redux_SS[i] = y_pan[i]
    

indices_SS=np.nonzero(x_redux_SS)[0]
#remove zeros
x_redux_SS=x_redux_SS[x_redux_SS!=0]
y_redux_SS=y_redux_SS[y_redux_SS!=0]  


O_isotope_redux_SS=del18Opan[indices_SS]#geochem.d18O[indices] 
Temp_redux_SS=temppan[indices_SS] #geochem.Temp[indices]
dip_angle=55 #dip of beds
A=[7659000,727000]; B=[7657000,729000] #First entry is northing, second easting
md=(A[0]-B[0])/(A[1]-B[1])#slope of the line
bd=A[0]-md*A[1] #y-intercept

#then, for every point, find its position (D) along that line, by defining a line P that is perpendicular
# to D
mp=-1/md #slope of this line is the same for all points
bp=np.subtract(y_redux_SS,mp*x_redux_SS) #y-intercept is unique for each point
#then, at each projection of point onto base line D, line D = line P
PosDx = (bd-bp)/(mp-md) #x-position in original axes
PosDy = mp*PosDx + bp #y-position in original axes
PosDtrue_SS = ((A[0]-PosDx)**2+(A[1]-PosDy)**2)**0.5#position along the line at the base of the formation, using A as the origin in the NW

#then, need to find the vertical distance above the line D for data points projected perpendicularly
tempx = np.subtract(PosDx,x_redux_SS)**2
tempy = np.subtract(PosDy,y_redux_SS)**2
DistOrig=np.add(tempx,tempy)**0.5 #Distance from point on line D to original point
 #calculate distance from line D to point Z, which is stratigraphic height
Z_SS=np.tan(np.deg2rad(dip_angle))*DistOrig


 #%%   
#--BK

x_redux_BK = np.zeros(len(x_pan))
y_redux_BK = np.zeros(len(y_pan)) 
for i in range(len(x_pan)):
    if x_pan[i]>728000 and x_pan[i]<733000 and y_pan[i]>7644000 and y_pan[i]<7649000:
        x_redux_BK[i] = x_pan[i]
        y_redux_BK[i] = y_pan[i]
    

indices_BK=np.nonzero(x_redux_BK)[0]
#remove zeros
x_redux_BK=x_redux_BK[x_redux_BK!=0]
y_redux_BK=y_redux_BK[y_redux_BK!=0]  


O_isotope_redux_BK=del18Opan[indices_BK]#geochem.d18O[indices] 
Temp_redux_BK=temppan[indices_BK] #geochem.Temp[indices]
dip_angle=55 #dip of beds
#A=[7645000,727000]; B=[7649000,731000] #First entry is northing, second easting
A=[763.988e4,722731]; B=[764.941e4,731385]
md=(A[0]-B[0])/(A[1]-B[1])#slope of the line
bd=A[0]-md*A[1] #y-intercept

#then, for every point, find its position (D) along that line, by defining a line P that is perpendicular
# to D
mp=-1/md #slope of this line is the same for all points
bp=np.subtract(y_redux_BK,mp*x_redux_BK) #y-intercept is unique for each point
#then, at each projection of point onto base line D, line D = line P
PosDx = (bd-bp)/(mp-md) #x-position in original axes
PosDy = mp*PosDx + bp #y-position in original axes
PosDtrue_BK = ((A[0]-PosDx)**2+(A[1]-PosDy)**2)**0.5#position along the line at the base of the formation, using A as the origin in the SW

#then, need to find the vertical distance above the line D for data points projected perpendicularly
tempx = np.subtract(PosDx,x_redux_BK)**2
tempy = np.subtract(PosDy,y_redux_BK)**2
DistOrig=np.add(tempx,tempy)**0.5 #Distance from point on line D to original point
 #calculate distance from line D to point Z, which is stratigraphic height
Z_BK=np.tan(np.deg2rad(dip_angle))*DistOrig

#--MW

x_redux_MW = np.zeros(len(x_pan))
y_redux_MW = np.zeros(len(y_pan)) 
for i in range(len(x_pan)):
    if x_pan[i]>725000 and x_pan[i]<735000 and y_pan[i]>7641500 and y_pan[i]<7646000:
        x_redux_MW[i] = x_pan[i]
        y_redux_MW[i] = y_pan[i]
    

indices_MW=np.nonzero(x_redux_MW)[0]
#remove zeros
x_redux_MW=x_redux_MW[x_redux_MW!=0]
y_redux_MW=y_redux_MW[y_redux_MW!=0]  


O_isotope_redux_MW=del18Opan[indices_MW]#geochem.d18O[indices] 
Temp_redux_MW=temppan[indices_MW] #geochem.Temp[indices]
dip_angle=55 #dip of beds
#A=[7642500,726500]; B=[7645500,727800] #First entry is northing, second easting
md=(A[0]-B[0])/(A[1]-B[1])#slope of the line
bd=A[0]-md*A[1] #y-intercept

#then, for every point, find its position (D) along that line, by defining a line P that is perpendicular
# to D
mp=-1/md #slope of this line is the same for all points
bp=np.subtract(y_redux_MW,mp*x_redux_MW) #y-intercept is unique for each point
#then, at each projection of point onto base line D, line D = line P
PosDx = (bd-bp)/(mp-md) #x-position in original axes
PosDy = mp*PosDx + bp #y-position in original axes
PosDtrue_MW = ((A[0]-PosDx)**2+(A[1]-PosDy)**2)**0.5#position along the line at the base of the formation, using A as the origin in the NW

#then, need to find the vertical distance above the line D for data points projected perpendicularly
tempx = np.subtract(PosDx,x_redux_MW)**2
tempy = np.subtract(PosDy,y_redux_MW)**2
DistOrig=np.add(tempx,tempy)**0.5 #Distance from point on line D to original point
 #calculate distance from line D to point Z, which is stratigraphic height
Z_MW=np.tan(np.deg2rad(dip_angle))*DistOrig

x_redux_A45 = np.zeros(len(x_pan))
y_redux_A45 = np.zeros(len(y_pan)) 
for i in range(len(x_pan)):
    if x_pan[i]>722000 and x_pan[i]<728000 and y_pan[i]>763000 and y_pan[i]<7643450:
        x_redux_A45[i] = x_pan[i]
        y_redux_A45[i] = y_pan[i]
    

indices_A45 = np.nonzero(x_redux_A45)[0]
#remove zeros
x_redux_A45 = x_redux_A45[x_redux_A45!=0]
y_redux_A45 = y_redux_A45[y_redux_A45!=0]  


O_isotope_redux_A45 = del18Opan[indices_A45]#geochem.d18O[indices] 
Temp_redux_A45=temppan[indices_A45] #geochem.Temp[indices]
dip_angle=55 #dip of beds
#A=[7642500,726500]; B=[7645500,727800] #First entry is northing, second easting
md=(A[0]-B[0])/(A[1]-B[1])#slope of the line
bd=A[0]-md*A[1] #y-intercept

#then, for every point, find its position (D) along that line, by defining a line P that is perpendicular
# to D
mp=-1/md #slope of this line is the same for all points
bp=np.subtract(y_redux_A45,mp*x_redux_A45) #y-intercept is unique for each point
#then, at each projection of point onto base line D, line D = line P
PosDx = (bd-bp)/(mp-md) #x-position in original axes
PosDy = mp*PosDx + bp #y-position in original axes
PosDtrue_A45 = ((A[0]-PosDx)**2+(A[1]-PosDy)**2)**0.5#position along the line at the base of the formation, using A as the origin in the NW

#then, need to find the vertical distance above the line D for data points projected perpendicularly
tempx = np.subtract(PosDx,x_redux_A45)**2
tempy = np.subtract(PosDy,y_redux_A45)**2
DistOrig=np.add(tempx,tempy)**0.5 #Distance from point on line D to original point
 #calculate distance from line D to point Z, which is stratigraphic height
Z_A45=np.tan(np.deg2rad(dip_angle))*DistOrig

#O_isotope_Whole = np.concatenate((O_isotope_redux_SS,O_isotope_redux,O_isotope_redux_BK,O_isotope_redux_MW,O_isotope_redux_A45))
O_isotope_Whole = np.concatenate((O_isotope_redux_SS,O_isotope_redux,O_isotope_redux_SE))

#Temp_Whole = np.concatenate((Temp_redux_SS,Temp_redux,Temp_redux_BK,Temp_redux_MW,Temp_redux_A45))
Temp_Whole = np.concatenate((Temp_redux_SS,Temp_redux,Temp_redux_SE))


#Z_Whole = np.concatenate((Z_SS,Z,Z_BK,Z_MW,Z_A45))
Z_Whole = np.concatenate((Z_SS,Z,Z_SE))

SS_true = np.subtract(PosDtrue_SS,PosDtrue_SS.min())
KC_true = np.add(np.subtract(PosDtrue,PosDtrue.min()),SS_true.max())
BK_true = np.add(np.subtract(PosDtrue_BK,PosDtrue_BK.min()),KC_true.max())
MW_true = np.add(np.subtract(PosDtrue_MW,PosDtrue_MW.min()),BK_true.max())
A45_true = np.add(np.subtract(PosDtrue_A45,PosDtrue_A45.min()),MW_true.max())
SE_true = np.add(np.subtract(PosDtrue_SE,PosDtrue_SE.min()),KC_true.max())
#X_Whole = np.concatenate((SS_true,KC_true,BK_true,MW_true,A45_true))
X_Whole = np.concatenate((SS_true,KC_true,SE_true))

Whole_water_initial = np.zeros(num_runspano)
Whole_water_outgoing = np.zeros(num_runspano)
Whole_moles_fluid = np.zeros(num_runspano)
Whole_W_R = np.zeros(num_runspano)
Whole_final_rock = np.zeros(num_runspano)
Whole_moles_Orock = np.zeros(num_runspano)
for irun in range(0,num_runspano):
    del18O_in_Whole = O_isotope_Whole
    Temp_Whole_in = Temp_Whole 
    X_Whole_in = X_Whole
    Z_Whole_in = Z_Whole
    
    num_takeout = np.random.randint(0,len(del18O_in_Whole))
#    num_replace = np.random.randint(0,len(del18O_in_MW)) 
#    del18O_in_MW[num_takeout] = del18O_in_MW[num_replace]
    
    del18O_in_Whole = np.delete(del18O_in_Whole,num_takeout)
    Temp_Whole_in =np.delete(Temp_Whole_in,num_takeout)
    X_Whole_in =  np.delete(X_Whole_in,num_takeout)
    Z_Whole_in =np.delete(Z_Whole_in,num_takeout)
   
    [Whole_iso_grid,Whole_x_grid,Whole_z_grid,conc_matrix_Whole,change_matrix_Whole,
     qsolved_Whole,moles_fluid_whole,W_R_taper,water_initial_Whole, water_outgoing_Whole, 
     temp_test_Whole,temp_grid_Whole,moles_Orock_Whole,Whole_epsilon_r_w] = \
    o_isotope_invert(del18O_in_Whole,Temp_Whole_in,X_Whole_in,Z_Whole_in,d18O_pano)
  
    Whole_water_initial[irun] = water_initial_Whole
    Whole_water_outgoing[irun] = water_outgoing_Whole
    Whole_moles_fluid[irun] = moles_fluid_whole
    Whole_W_R[irun] = W_R_taper
    Whole_final_rock[irun] = del18O_in_Whole.mean()
    Whole_moles_Orock[irun] = moles_Orock_Whole
    
#%%Hydration sensitivity test    

#using relationship from Muehlenbachs and Clayton, 1972, where dRf = 0.97*wt% water + dRi
# so use gridded final isotope values to predict water %. Then, convert wt.% to 
#moles O in hydrated minerals, divide by moles O in the rock to get a "hydration" F/R ratio
# this should be a lot smaller than what the inversion predicts for total F/R, so 
#just subtract hydration F/R from inversion F/R, calculate new estimated dWi to test
# how much accounting for hydration affects this value
# For grid points that are depleted comared to starting dRi (5.8 permil), we
# assume an average water content equal to the average of all enriched rock value locations

inverse_water_initial = [EPR_water_initial.mean(),ohmoto_water_initial.mean(),solea_water_initial.mean(),Whole_water_initial.mean()]    
iso_grids_all = [EPR_iso_grid,ohmoto_iso_grid,solea_iso_grid,Whole_iso_grid]    
final_rock_all = [EPR_final_rock.mean(),ohmoto_final_rock.mean(),solea_final_rock.mean(),Whole_final_rock.mean()]
rock_O_all = [moles_Orock_EPR,moles_Orock_ohmoto,moles_Orock_solea,moles_Orock_Whole]
W_R_all = [EPR_W_R.mean(),ohmoto_W_R.mean(),solea_W_R.mean(),Whole_W_R.mean()]
Delt_RW = [np.mean(EPR_epsilon_r_w),np.mean(ohmoto_epsilon_r_w),np.mean(solea_epsilon_r_w),np.mean(Whole_epsilon_r_w)]
hydro_W_R = np.zeros(len(inverse_water_initial))
hydro_water_initial = np.zeros(len(inverse_water_initial))
test_initial = np.zeros(len(inverse_water_initial))
dRi = [d18Oinit_morb,d18Oinit_hoku,d18Oinit_morb,d18Oinit_morb]
for iloc in range(0,len(inverse_water_initial)):
    m = 0.97; 
    water = (iso_grids_all[iloc] - dRi[iloc])/m
#    np.where(wt_per_EPR<0, 0.1, wt_per_EPR)
    enriched = np.where(water>0); depleted = np.where(water<0)
    enriched_mean = np.mean(water[enriched])
    water[depleted] = enriched_mean; water = (water/100)
    moles_O_hydrated = water.mean()*3/(18*0.086)*rock_O_all[iloc]
    hydro_W_R[iloc] = moles_O_hydrated/rock_O_all[iloc]
    aa = final_rock_all[iloc] - dRi[iloc]; bb = W_R_all[iloc]-hydro_W_R[iloc]
    hydro_water_initial[iloc] = aa/bb + final_rock_all[iloc] - Delt_RW[iloc]
    test_initial[iloc] = aa/W_R_all[iloc] + final_rock_all[iloc] - Delt_RW[iloc]
    
#    (d18Of-d18Oinit)/(W_R_predicted) + d18Of - np.nanmean(epsilon_r_w)
#%% Seawater oxygen isotope exchange model
t = 4.4 #time in Gyr

Wo = 7 #original seawater d18O
W_SS = -1 #steady state 
num_steps = 100 #num of initial model steps
time = np.linspace(0,4.5,num=num_steps) #sample every 250 myr
weath_time_on = 4.5-2.8#in Ga
weath_time_early = 4.5-4.43
weath_time_late = 4.5-0.9
# rate constants in Gyr-1, from Muehlenbachs, 1998
k_weath = 8 #nominal 8continental weathering
k_growth = 1.2 #nominal 1.2continental growth
k_hiT = 14.6 #nominal 14.6high temperature seafloor 
k_loT = 1.7 #nominal 1.7low temp seafloor/seafloor weathering 
k_W_recycling = 0.6 #nominal 0.6water recycling at subduction zones

#fractionations (permil) btwn rock and water, from Muehlenbachs, 1998 except weathering, which we tuned to reproduce -1permil ocean

Delt_weath =  13  #mueh = 9.6 nominal 13, newer 17
Delt_growth =  9.8 #mueh = 9.8 
Delt_hiT_mid =  1.5 #meuh = 4.1
Delt_hiT = 4.1 
#Delt_hiT_mid = np.zeros(len(time))

Delt_lowT =  9.3  #mueh = 9.3 
Delt_water_recycling = 2.5 # mueh = 2.5 
    

#calculate steady state in 250 myr increments 
del_graniteo = np.linspace(7.8,7.8,num=num_steps)

del_basalto = 5.5
del_WR = 7
bb= 0.1


Delt_hiT_change = (Delt_hiT-Delt_hiT_mid)+(Delt_hiT-Delt_hiT_mid)*0.5*(1+np.tanh((np.subtract(time,weath_time_on)/bb))) -1
Delt_hiT_change_late = (Delt_hiT-Delt_hiT_mid)+(Delt_hiT-Delt_hiT_mid)*0.5*(1+np.tanh((np.subtract(time,weath_time_late)/bb))) -1
k_growth_change = 0.5*k_growth*(1+np.tanh((np.subtract(time,weath_time_on)/bb)))
k_growth_late = 0.5*k_growth*(1+np.tanh((np.subtract(time,weath_time_late)/bb)))
k_growth_early = 0.5*k_growth*(1+np.tanh((np.subtract(time,weath_time_early)/bb)))

k_weathering_change =0.5*k_weath*(1+np.tanh((np.subtract(time,weath_time_on)/bb)))
k_weathering_late =0.5*k_weath*(1+np.tanh((np.subtract(time,weath_time_late)/bb)))
k_weathering_early =0.5*k_weath*(1+np.tanh((np.subtract(time,weath_time_early)/bb)))

k_loT_change = k_loT*np.ones(time.size) #keep it the same
k_hiT_change = k_hiT*np.ones(time.size)

k_water_change = k_W_recycling*np.ones(time.size) #


del_steady_change = np.zeros(time.size)
del_steady_early = np.zeros(time.size)
del_steady_late = np.zeros(time.size)
k_sum = np.zeros(time.size)
k_sum_early = np.zeros(time.size)
k_sum_late = np.zeros(time.size)

for istep in range(0,time.size):
    top = np.sum([k_weathering_change[istep]*(del_graniteo[istep]-Delt_weath),\
                  k_growth_change[istep]*(del_graniteo[istep]-Delt_growth),\
                  k_hiT_change[istep]*(del_basalto-Delt_hiT_change[istep]),\
                  k_loT_change[istep]*(del_basalto-Delt_lowT),\
                  k_water_change[istep]*(del_WR-Delt_water_recycling)])
    top_early = np.sum([k_weathering_early[istep]*(del_graniteo[istep]-Delt_weath),\
                  k_growth_early[istep]*(del_graniteo[istep]-Delt_growth),\
                  k_hiT_change[istep]*(del_basalto-Delt_hiT),\
                  k_loT_change[istep]*(del_basalto-Delt_lowT),\
                  k_water_change[istep]*(del_WR-Delt_water_recycling)])
    top_late = np.sum([k_weathering_late[istep]*(del_graniteo[istep]-Delt_weath),\
                  k_growth_late[istep]*(del_graniteo[istep]-Delt_growth),\
                  k_hiT_change[istep]*(del_basalto-Delt_hiT_change_late[istep]),\
                  k_loT_change[istep]*(del_basalto-Delt_lowT),\
                  k_water_change[istep]*(del_WR-Delt_water_recycling)])
    k_sum[istep] = np.sum([k_weathering_change[istep],k_growth_change[istep],k_hiT_change[istep],k_loT_change[istep],k_water_change[istep]])
    k_sum_early[istep] = np.sum([k_weathering_early[istep],k_growth_early[istep],k_hiT_change[istep],k_loT_change[istep],k_water_change[istep]])
    k_sum_late[istep] = np.sum([k_weathering_late[istep],k_growth_late[istep],k_hiT_change[istep],k_loT_change[istep],k_water_change[istep]])
    
    del_steady_change[istep] = top/k_sum[istep]
    del_steady_early[istep] = top_early/k_sum_early[istep]
    del_steady_late[istep] = top_late/k_sum_late[istep]
       
#calculate dW at for each steady state
time_new = np.linspace(0.01, 4.5,num=1e3)
f1 = sp.interpolate.interp1d(time,del_steady_change)
f2 = sp.interpolate.interp1d(time,k_sum)
steady_interp = f1(time_new)
k_sum_interp = f2(time_new)
f1_late = sp.interpolate.interp1d(time,del_steady_late)
f2_late = sp.interpolate.interp1d(time,k_sum_late)
steady_interp_late = f1_late(time_new)
k_sum_interp_late = f2_late(time_new)

steady_interp_late = f1_late(time_new)
k_sum_interp_late = f2_late(time_new)

f1_early = sp.interpolate.interp1d(time,del_steady_early)
f2_early = sp.interpolate.interp1d(time,k_sum_early)
steady_interp_early = f1_early(time_new)
k_sum_interp_early = f2_early(time_new)

dW_middle = np.add(np.subtract(Wo,steady_interp)*np.exp(-np.multiply(time_new,k_sum_interp)),steady_interp) 
dW_early = np.add(np.subtract(Wo,steady_interp_early)*np.exp(-np.multiply(time_new,k_sum_interp_early)),steady_interp_early) 
dW_late = np.add(np.subtract(Wo,steady_interp_late)*np.exp(-np.multiply(time_new,k_sum_interp_late)),steady_interp_late) 

decay_const_low = 0.02
decay_const_high = 0.04
dW_decay_low = []#np.zeros(len(time_new))
dW_decay_high = []#np.zeros(len(time_new))
whatstep=[]



for istep in range(0,len(time_new)):
    if time_new[istep]<=1.5:
        dW_decay_low.append(np.add(np.subtract(Wo,steady_interp[-1])*np.exp(-np.multiply(time_new[istep],decay_const_low*k_sum[-1])),steady_interp[-1]))
        dW_decay_high.append(np.add(np.subtract(Wo,steady_interp[-1])*np.exp(-np.multiply(time_new[istep],decay_const_high*k_sum[-1])),steady_interp[-1]))
        temp_Wo_low = dW_decay_low[istep]
        temp_Wo_high = dW_decay_high[istep]
        whatstep.append(istep)
knickpoint = whatstep[-1]
time_late = np.flip(time_new[-1] - time_new[knickpoint+1:],0) 
low_test = []
high_test = []
for istep in range(0,len(time_late)):
    dW_decay_low.append(np.add(np.subtract(temp_Wo_low,steady_interp[-1])*np.exp(-np.multiply(time_late[istep],0.4*k_sum[-1])),steady_interp[-1]))
    dW_decay_high.append(np.add(np.subtract(temp_Wo_high,steady_interp[-1])*np.exp(-np.multiply(time_late[istep],0.5*k_sum[-1])),steady_interp[-1]))
    low_test.append(np.add(np.subtract(temp_Wo_low,steady_interp[-1])*np.exp(-np.multiply(time_late[istep],0.03*k_sum[-1])),steady_interp[-1]))



#%%
#load in Lear et al., 2000 ocean d18O estimate
lear_data=pd.read_csv('lear_ocean_d18O.csv')    
lear_d18O = pd.Series(lear_data['d18O']).values
lear_age = pd.Series(lear_data['age_ma']).values 
lear_age_Ga = np.subtract(4.5,np.divide(lear_age,1e3))

#%% Save output to LaTeX table
#calculate Wing and Ferry equilibrium stats without outliers
inverse_n =[wing_water_initial.size,hoku50_water_initial.size,skaer_water_initial.size,EPR_water_initial.size,
                ohmoto_water_initial.size,solea_water_initial.size,Whole_water_initial.size]
wing_water_filter = wing_water_initial[wing_water_initial < 4]

temp_avgs = np.rint([temp_wing.mean(),temp_hoku50.mean(),temp_skaer.mean(),temp.mean(),
                    temp_ohmoto.mean(),temp_solea.mean(),temppan.mean()])
 
temp_avgs=temp_avgs.astype(int)

    
inverse_avgs = np.around([wing_water_filter.mean(),hoku50_water_initial.mean(),skaer_water_initial.mean(),EPR_water_initial.mean(),
                ohmoto_water_initial.mean(),solea_water_initial.mean(),Whole_water_initial.mean()], decimals=2)
inverse_std =    np.around([wing_water_filter.std(),hoku50_water_initial.std(),skaer_water_initial.std(),EPR_water_initial.std(),
                ohmoto_water_initial.std(),solea_water_initial.std(),Whole_water_initial.std()], decimals=2)
inverse_table = ['a']*inverse_avgs.size

inverse_avgs_str = [str(i) for i in (inverse_avgs.tolist())]
inverse_std_str = [str(i) for i in (inverse_std.tolist())]
pm='\pm'
for i in range(len(inverse_avgs)):
    inverse_table[i] = '$'+inverse_avgs_str[i]+pm+inverse_std_str[i]+'$'

moles_fluid_avgs = [wing_moles_fluid.mean(),hoku50_moles_fluid.mean(),skaer_moles_fluid.mean(),EPR_moles_fluid.mean(),
                    ohmoto_moles_fluid.mean(),solea_moles_fluid.mean(),Whole_moles_fluid.mean()]
moles_fluid_std = [wing_moles_fluid.std(),hoku50_moles_fluid.std(),skaer_moles_fluid.std(),EPR_moles_fluid.std(),
                    ohmoto_moles_fluid.std(),solea_moles_fluid.std(),Whole_moles_fluid.std()]


import math
def orderOfMagnitude(number):
    return math.floor(math.log(number, 10))

fluid_table = ['a']*len(moles_fluid_avgs)
for i in range(len(moles_fluid_avgs)):
    oom = orderOfMagnitude(moles_fluid_avgs[i])
    aa = np.round(moles_fluid_avgs[i]/10**(oom),decimals=1)
    aa=str(aa)
    bb = np.round(moles_fluid_std[i]/10**(oom),decimals=2)
    bb=str(bb)
    oom=str(oom)
    power = '\times10^{'+oom+'}'
    fluid_table[i] = '$'+aa+pm+bb+power+'$'
    
moles_rock_avgs = [wing_moles_Orock.mean(),hoku50_moles_Orock.mean(),skaer_moles_Orock.mean(),EPR_moles_Orock.mean(),
                    ohmoto_moles_Orock.mean(),solea_moles_Orock.mean(),Whole_moles_Orock.mean()]
moles_rock_std = [wing_moles_Orock.std(),hoku50_moles_Orock.std(),skaer_moles_Orock.std(),EPR_moles_Orock.std(),
                    ohmoto_moles_Orock.std(),solea_moles_Orock.std(),Whole_moles_Orock.std()]

Orock_table = ['a']*len(moles_rock_avgs)
for i in range(len(moles_rock_avgs)):
    oom = orderOfMagnitude(moles_rock_avgs[i])
    aa = np.round(moles_rock_avgs[i]/10**(oom),decimals=1)
    aa=str(aa)
    bb = np.round(moles_rock_std[i]/10**(oom),decimals=2)
    bb=str(bb)
    oom=str(oom)
    power = '\times10^{'+oom+'}'
    Orock_table[i] = '$'+aa+pm+bb+power+'$'    
    
  

W_R_avgs = np.around([wing_W_R.mean(),hoku50_W_R.mean(),skaer_W_R.mean(),EPR_W_R.mean(),
                ohmoto_W_R.mean(),solea_W_R.mean(),Whole_W_R.mean()], decimals=2)
W_R_std =    np.around([wing_W_R.std(),hoku50_W_R.std(),skaer_W_R.std(),EPR_W_R.std(),
                ohmoto_W_R.std(),solea_W_R.std(),Whole_W_R.std()], decimals=2)
W_R_table = ['a']*W_R_avgs.size

W_R_avgs_str = [str(i) for i in (W_R_avgs.tolist())]
W_R_std_str = [str(i) for i in (W_R_std.tolist())]
pm='\pm'
for i in range(len(W_R_avgs)):
    W_R_table[i] = '$'+W_R_avgs_str[i]+pm+W_R_std_str[i]+'$'

AA = -3.9e3; BB = 3.42e6 #1000ln alpha(rock-water) = AA/T + BB/T^2 
temp_kelvin = np.add(temp_avgs,273)

Big_delt = np.around([np.mean(wing_epsilon_r_w),np.mean(hoku50_epsilon_r_w),np.mean(skaer_epsilon_r_w),np.mean(EPR_epsilon_r_w),
                np.mean(ohmoto_epsilon_r_w),np.mean(solea_epsilon_r_w),np.mean(Whole_epsilon_r_w)], decimals=2)

Rock_initials = [d18Oinit_wing,d18Oinit_hoku50,d18Oinit_skaer,d18Oinit_morb,d18Oinit_ohmoto,d18Oinit_solea,d18O_pano]
Rock_finals = np.around([wing_final_rock.mean(),hoku50_final_rock.mean(),skaer_final_rock.mean(),EPR_final_rock.mean(),
                ohmoto_final_rock.mean(),solea_final_rock.mean(),Whole_final_rock.mean()], decimals=1)

ocean_del18O = ['-','-',]    
column_names = ['\textbf{Setting$^{\text{ref}}$}', '\textbf{$n$}','\textbf{Temp ($^\circ$C)}',
                                 '\textbf{Incoming $\delta^{18}$O}','\textbf{Inverse $\delta^{18}$O}',
                                 '\textbf{Moles O Fluid}','\textbf{Moles O Rock}','\textbf{F/R}','\textbf{$\Delta$}','\textbf{$\delta R_i$}',
                                 '\textbf{$\delta R_f$}']

inverse_output = {column_names[0]: ['Equilibrium$^1$','Kinetic$^2$','Skeargaard$^3$','East Pacific Rise$^4$',
                  'Fukazawa$^5$','Solea$^6$','Panorama$^7$'],
        column_names[1]: inverse_n,
        column_names[2]:temp_avgs,
        column_names[3]: [0,0,-14,0.2,-0.55,-1,'-'],
        column_names[8]: Big_delt,
        column_names[9]: Rock_initials,
        column_names[10]: Rock_finals,
        column_names[5]: fluid_table,
        column_names[6]: Orock_table,
        column_names[7]: W_R_table,
        column_names[4]: inverse_table
        }

df = pd.DataFrame(inverse_output,columns= column_names)

bot_refs = '\multicolumn{5}{c}{1-\cite{Wing_and_Ferry_2007},2-\cite{Cathles_1983},3-\cite{Norton_and_Taylor_1979},4-\cite{Gillis_et_al_2001},5-\cite{Green_et_al_1983},6-\cite{Schiffman_and_Smith_1988},7-\cite{Brauhart_et_al_2000}}'
df_tex = df.to_latex(index=False,escape=False,column_format='l c c c c c c c c c c')
df_tex=df_tex.replace('toprule','hline').replace('bottomrule','hline').replace('midrule','hline')

#add multicolumn row with refs after bottom hline
h_idx = df_tex.rfind('\hline')+6
df_tex = df_tex[:h_idx]+bot_refs+df_tex[h_idx:]


file = open('inverse_output.tex','w')  
file.write(df_tex) 
file.close() 

model_column_names = ['\textbf{Flux}', '\multicolumn{1}{c}{\textbf{$k_i$} (Gyr$^{-1}$)}','\multicolumn{1}{c}{\textbf{$\Delta$} ($\permil$) } ',
                                 '\multicolumn{1}{c}{\textbf{$\delta^o$} ($\permil$) }']
hiT_cols = [np.round(Delt_hiT_change.min(),decimals=1),np.round(Delt_hiT_change.max(),decimals=1)]
hiT_cols = np.str(hiT_cols[0])+ ' - ' +np.str(hiT_cols[1])
ks=np.round([k_weath,k_growth,k_hiT,k_loT,k_W_recycling],decimals=1)
model_conditions = {model_column_names[0]: ['Continental Weathering','Continental Recycling','High Temperature alteration',
                  'Low Temperature Alteration','Water Recycling'],
        model_column_names[1]: ks,
        model_column_names[2]:[Delt_weath,Delt_growth,hiT_cols,Delt_lowT,Delt_water_recycling],
        model_column_names[3]: [del_graniteo[0],del_graniteo[0],del_basalto,del_basalto,del_WR],
        }

df_model = pd.DataFrame(model_conditions,columns= model_column_names)

df_model_tex = df_model.to_latex(index=False,escape=False,column_format='l . . . ')

df_model_tex=df_model_tex.replace('toprule','hline').replace('bottomrule','hline').replace('midrule','hline')
file = open('d18O_conditions.tex','w')  
file.write(df_model_tex) 
file.close() 

#hydration sensitivity table
hydr_column_names= [ '\textbf{Setting}','\multicolumn{1}{c}{\textbf{$\delta R_i$}}','\multicolumn{1}{c}{\textbf{F/R}}','\multicolumn{1}{c}{\textbf{Hydrated F/R}}',
                    '\multicolumn{1}{c}{\textbf{Hydrated $\delta^{18}$O}}','\multicolumn{1}{c}{\textbf{Inverse $\delta^{18}$O}}',
                    '\multicolumn{1}{c}{\textbf{Difference}}']
settings = ['EPR','Fukazawa','Solea','Panorama']
hydro_wat_table = np.around(hydro_water_initial,decimals=2)
inverse_init_table = np.round(inverse_water_initial,decimals=2)
difference = np.round([inverse_water_initial - hydro_wat_table],decimals=2)
difference = np.reshape(difference,[4,])
hydro_W_R_table = np.round(hydro_W_R,decimals=3)
W_R_table[-4:-1]

hydro_sensitivity = {hydr_column_names[0]: ['East Pacific Rise','Fukazawa','Solea','Panorama'],
        hydr_column_names[1]: Rock_initials[-4:],
        hydr_column_names[2]: W_R_table[-4:],
        hydr_column_names[3]:hydro_W_R_table,
        hydr_column_names[4]: hydro_wat_table,
        hydr_column_names[5]: inverse_init_table,
        hydr_column_names[6]: difference
        }
df_sens = pd.DataFrame(hydro_sensitivity,columns= hydr_column_names)

df_sens_tex = df_sens.to_latex(index=False,escape=False,column_format='l c c c c c c')

df_sens_tex=df_sens_tex.replace('toprule','hline').replace('bottomrule','hline').replace('midrule','hline')
file = open('hydro_sensitivity.tex','w')  
file.write(df_sens_tex) 
file.close()
#%% Figures 
plt.close('all')

isogrid_cmap = plt.cm.BrBG
tempgrid_cmap = plt.cm.coolwarm
fig1=plt.figure(num=1); plt.clf(); #plt.title('grid')
plt.subplot(1,2,1)
im = plt.imread('geo_map2.png')
implot = plt.imshow(im)

map_lims = np.divide([718500,735000, 7637000, 7661000],1e4)
plt.imshow(im, zorder=0, extent=map_lims)

plt.xlim([map_lims[0],map_lims[1]]); plt.ylim([map_lims[2],map_lims[3]])
plt.plot(x_pan/1e4,y_pan/1e4,'wo',markeredgecolor='k',markersize=7)
plt.text(72.3,765.85,'A')
plt.text(72.3,763.8,"A'") 
plt.legend(['Sample Sites'],loc='upper left',bbox_to_anchor=(0,1.1))
ax=plt.gca()
#xlabels = np.divide(plt.xticks()[0],1e4).tolist()
labels = [item.get_text() for item in ax.get_xticklabels()]
#labels = xlabels

#ax.set_xticklabels(xlabels)
#ax.set_yticklabels(empty_string_labels)
plt.xlabel('Easting ($10^4$ m)'); plt.ylabel('Northing ($10^4$ m)')

plt.subplot(2,2,2)
vert_ticks = np.linspace(0,3000,num=5)
plt.ylim([0,Z_Whole.max()])
plt.yticks(vert_ticks)
#plt.contourf(Whole_iso_grid,cmap=isogrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()]) #Whole_x_grid,Whole_z_grid,
numsteps=50
xi = np.arange(X_Whole.min(),X_Whole.max(),(X_Whole.max()-X_Whole.min())/numsteps)
yi = np.arange(Z_Whole.min(),Z_Whole.max(),(Z_Whole.max()-Z_Whole.min())/numsteps)
xi,yi = np.meshgrid(xi,yi)  
from scipy.interpolate import griddata 
zi = griddata((X_Whole,Z_Whole),O_isotope_Whole,(xi,yi),method='linear')
levels = np.linspace(2,15,10)
#plt.contourf(O_isotope_grid,cmap=isogrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
plt.contourf(xi,yi,zi,cmap=isogrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])

plt.plot(X_Whole,Z_Whole,'w.',markeredgecolor='k');
plt.plot(Whole_x_grid,Whole_z_grid,'k.')
plt.ylabel('Vertical (m)')
ax1=plt.gca()
#ax1.set_aspect('equal')
#ax1.set_xticklabels('')
ax1.xaxis.tick_top()
#plt.plot(x_redux_A45,y_redux_A45,'w^',markeredgecolor='k',markersize=12)
#cm = plt.cm.get_cmap('RdYlBu')
cbar=plt.colorbar(orientation='horizontal',pad=0.05)
#cbar.ax.tick_params(axis='x',direction='in',labeltop='on')
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_xlabel('$\delta^{18}$O (‰ VSMOW)')
iso_ticks = [3,5,7,9,11,13,15]
cbar.set_ticks(iso_ticks)
cbar.set_ticklabels(iso_ticks)

plt.subplot(2,2,4)

#plt.contourf(temp_grid_Whole,cmap=tempgrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
zii = griddata((X_Whole,Z_Whole),Temp_Whole,(xi,yi),method='linear')
levels = np.linspace(80,360,7)
plt.contourf(xi,yi,zii,levels,cmap=tempgrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
plt.plot(X_Whole,Z_Whole,'w.',markeredgecolor='k')
plt.plot(Whole_x_grid,Whole_z_grid,'k.')
plt.ylim([0,Z_Whole.max()])
plt.yticks(vert_ticks)
plt.ylabel('Vertical (m)')
plt.xlabel('Horizontal (m)');
plt.text(X_Whole.min(),Z_Whole.min()-450,'A')
plt.text(0.95*X_Whole.max(),Z_Whole.min()-450,"A'") 
ax2 = plt.gca()
ax2.xaxis.tick_top()
ax2.set_xticklabels('')
#ax2.set_yticklabels('')
#ax2.set_aspect('equal')
#plt.plot(x_redux_A45,y_redux_A45,'w^',markeredgecolor='k',markersize=12)
cbar=plt.colorbar(orientation='horizontal',pad=0.2)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_xlabel('Temp. ($^\circ$C)')

plt.show()
fig1.set_size_inches(14, 8.5)    
fig1.savefig('iso_temp_contours.eps',metadata='eps')


num_bins = 10

alphabar = 0.8
EPR_color=np.divide([51,34,136],255)
ohmoto_color=np.divide([221,204,119],255) 
solea_color=np.divide([204,102,119],255)
pano_color=np.divide([136,204,238],255)
wing_color=np.divide([68,170,153],255)
hoku50_color=np.divide([170,68,153],255)
skaer_color=np.divide([17,119,51],255)

EPR_weights = np.ones_like(EPR_water_initial)/float(len(EPR_water_initial))
ohmoto_weights = np.ones_like(ohmoto_water_initial)/float(len(ohmoto_water_initial))
solea_weights = np.ones_like(solea_water_initial)/float(len(solea_water_initial))
pano_weights = np.ones_like(Whole_water_initial)/float(len(Whole_water_initial))

wing_weights = np.ones_like(wing_water_initial)/float(len(wing_water_initial))
hoku50_weights = np.ones_like(hoku50_water_initial)/float(len(hoku50_water_initial))
skaer_weights = np.ones_like(skaer_water_initial)/float(len(skaer_water_initial))

#plt.hist(EPR_water_initial, weights=weights)

fig1s, ax1 = plt.subplots()  

#plt.subplot(1,2,1)
binwidth = .8
ax1.hist(wing_water_initial,bins=np.arange(min(wing_water_initial), max(wing_water_initial) + binwidth, binwidth), weights=wing_weights,color=wing_color,ec='k',alpha=0.5)
ax1.hist(hoku50_water_initial,bins=np.arange(min(hoku50_water_initial), max(hoku50_water_initial) + binwidth, binwidth), weights=hoku50_weights,color=hoku50_color,ec='k',alpha=0.5)
ax1.hist(skaer_water_initial,bins=np.arange(min(skaer_water_initial), max(skaer_water_initial) + binwidth, binwidth), weights=skaer_weights,color=skaer_color,ec='k',alpha=0.5)
#ax1=plt.gca()
ax1.set_title('Forward models')
ax1.set_xlabel('$\delta^{18}$O incoming (‰ VSMOW)'); ax1.set_ylabel('Fraction of runs')

ax1.annotate("",
            xy=(-14, 0.78), xycoords='data',
            xytext=(-14, 0.93), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
ax1.annotate("",
            xy=(0, 1), xycoords='data',
            xytext=(0, 1.05), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )            

ax1.legend(['Equilibrium','Forward','Skaergaard'])

axins = inset_axes(ax1,
                   width="30%",  # width = 30% of parent_bbox
                   height=1.3,  # height : 1 inch
                   loc='center')
binwidth = .1
axins.hist(wing_water_initial,bins=np.arange(min(wing_water_initial), max(wing_water_initial) + binwidth, binwidth), weights=wing_weights,color=wing_color,ec='k',alpha=0.5)
axins.hist(hoku50_water_initial,bins=np.arange(min(hoku50_water_initial), max(hoku50_water_initial) + binwidth, binwidth), weights=hoku50_weights,color=hoku50_color,ec='k',alpha=0.5)
axins.set_xlim(-1,1)
axins.set_ylabel('Fraction of runs'); axins.set_xlabel('$\delta^{18}$O incoming (‰ VSMOW)')

fig1s.set_size_inches(12, 8.5)    
fig1s.savefig('forward_results.eps',metadata='eps')

fig2, ax2 =plt.subplots()
binwidth = .1
ax2 = plt.gca()
ax2.hist(EPR_water_initial,bins=np.arange(min(EPR_water_initial), max(EPR_water_initial) + binwidth, binwidth), weights=EPR_weights,color=EPR_color,ec='k',alpha=alphabar)
ax2.hist(ohmoto_water_initial,bins=np.arange(min(ohmoto_water_initial), max(ohmoto_water_initial) + binwidth, binwidth), weights=ohmoto_weights,color=ohmoto_color,ec='k',alpha=alphabar)
ax2.hist(solea_water_initial,bins=np.arange(min(solea_water_initial), max(solea_water_initial) + binwidth, binwidth), weights=solea_weights,color=solea_color,ec='k',alpha=alphabar)
ax2.hist(Whole_water_initial,bins=np.arange(min(Whole_water_initial), max(Whole_water_initial) + binwidth, binwidth), weights=pano_weights,color=pano_color,ec='k',alpha=alphabar)

#ax2.set_title('Crustal Sections')
#
#ax2.set_yticklabels('')
ax2.legend(['EPR (0.2 Ma)','Fukazawa (14.2 Ma)','Solea (90 Ma)','Panorama (3.24 Ga)'],loc=[0.4,0.7]) #,bbox_to_anchor=(0.55, 1)
ax2.set_xlabel('$\delta^{18}$O incoming (‰ VSMOW)'); ax2.set_ylabel('Fraction of runs')


fig2.set_size_inches(12, 8.5)    
fig2.savefig('inverse_results.eps',metadata='eps')
 
#plt.figure(num=2)
#binwidth = .1
#plt.hist(wing_water_initial,bins=np.arange(min(wing_water_initial), max(wing_water_initial) + binwidth, binwidth), weights=wing_weights,color=wing_color,ec='k',alpha=0.5)
#plt.hist(hoku50_water_initial,bins=np.arange(min(hoku50_water_initial), max(hoku50_water_initial) + binwidth, binwidth), weights=hoku50_weights,color=hoku50_color,ec='k',alpha=0.5)
#plt.xlim([-1.5,1.5])
##plt.legend(['Wing and Ferry','Cathles'])

#plt.xlabel('$\delta^{18}$O incoming'); plt.ylabel('Fraction of runs')



#%%
#Ocean O-isotope
inversion_ages = [4.45,3.2,0.0916,0.01432,0.002]
inv_age_err = [0,0.005,0.0014,0.0005,0]
inversion_d18O = [6.5,Whole_water_initial.mean(),solea_water_initial.mean(),\
                  ohmoto_water_initial.mean(),EPR_water_initial.mean()]
inversion_err = np.multiply(2,[0.25,Whole_water_initial.std(),solea_water_initial.std(),\
                  ohmoto_water_initial.std(),EPR_water_initial.std()])

previous_ages = [3.8,0.76,0.443,0.065] #Pope, Hodel, Muehlenbachs, Gregory and Taylor
previous_d18O = [2.3,-1.33,0,-0.4]
previous_err = [1.7,.98,2,1]
time_labels = ['4.5','4','3.5','3','2.5','2','1.5','1','0.5','0']
time_ticks = np.linspace(0,4.5,10)

#fig3, (ax1) = plt.subplots(1, 1, sharey=True)
#ax1 = plt.subplot(1,1,1)
#ax1.plot(time,k_growth_change,'k:')
#ax1.plot(time,k_weathering_change,'k--')
#ax1.plot(time,k_loT_change,'k-.')
#ax1.plot(time,k_hiT_change,'k-')
#ax1.plot(time,k_water_change,'k.-')
#ax1.legend(['Continental recycling','Continental weathering','low T alteration','high T alteration','Water recycling']\
#           ,bbox_to_anchor=(0.01, 0.93), loc=2, borderaxespad=0.,fontsize=8)
#plt.xlabel('Age (Ga)'); plt.ylabel('Seawater $\delta^{18}$O')
#locs, labels = plt.xticks()           # Get locations and labels
#
#plt.xticks(time_ticks, time_labels)  # Set locations and labels
#plt.xticks(time_ticks, time_labels)  
#ax1.set_xlim([0,4.5])
#plt.ylabel('Rate (Gyr $^{-1}$)')
this_study_color='xkcd:melon'
#fig3, (ax2) = plt.subplots(1, 1, sharey=True)
fig3 = plt.figure();plt.clf()
#ax2 = plt.subplot(1,1,1) # subplot(2,1,2) is now active
#ax2.set_xlim([0,4.5])
plt.xticks(time_ticks, time_labels)  # Set locations and labels
plt.xlabel('Age (Ga)'); plt.ylabel('Seawater $\delta^{18}$O (‰ VSMOW)')
#plt.yticks([])
#plt.title('Seawater $\delta^{18}$O constraints')

plt.plot(time_new,dW_decay_low,' ',linewidth=2, label=str())
plt.plot(time_new,dW_decay_high,' ',linewidth=2, label=str())
lgd1=plt.fill_between(time_new, dW_decay_low, dW_decay_high, color='xkcd:teal',alpha=0.5)


#ax2.plot(time_new,dW_late,'k--',linewidth=2)
lgd2=plt.errorbar(np.subtract(4.5,previous_ages),previous_d18O,yerr=previous_err,fmt='o',marker='o', mfc='xkcd:silver',\
         mec='xkcd:silver',ecolor='xkcd:silver',capsize=10, elinewidth=20)

lgd3=plt.errorbar(np.subtract(4.5,inversion_ages[1:]),inversion_d18O[1:],yerr=inversion_err[1:],fmt='o',marker='o', mfc=this_study_color,\
         mec=this_study_color,ecolor=this_study_color,capsize=10, elinewidth=2)

plt.errorbar(np.subtract(4.5,inversion_ages[0]),inversion_d18O[0],yerr=inversion_err[0],fmt='o',marker='o', mfc=this_study_color,\
         mec=this_study_color,ecolor=this_study_color,capsize=10, elinewidth=20)

lgd4,=plt.plot(time_new,dW_early,'k-.',linewidth=2)
lgd5,=plt.plot(time_new,dW_middle,'k',linewidth=2)

import matplotlib.patches as patches
#a_em =patches.Rectangle((4.5-2.5,-3),.03,10,linewidth=1,edgecolor='k',facecolor='xkcd:cloudy blue',alpha=0.4) 
#ax2.add_patch(a_em)
#b_em =patches.Rectangle((4.5-0.7,-3),.13,10,linewidth=1,edgecolor='k',facecolor='xkcd:cloudy blue',alpha=0.4) 
#ax2.add_patch(b_em)
#c_em =patches.Rectangle((4.5-3.5,-3),.03,10,linewidth=1,edgecolor='k',facecolor='xkcd:cloudy blue',alpha=0.4) 
#ax2.add_patch(c_em)
#d_em =patches.Rectangle((4.5-3.2,-3),.03,10,linewidth=1,edgecolor='k',facecolor='xkcd:cloudy blue',alpha=0.4) 
#ax2.add_patch(d_em)
#e_em =patches.Rectangle((4.5-3,-3),.03,10,linewidth=1,edgecolor='k',facecolor='xkcd:cloudy blue',alpha=0.4) 
#ax2.add_patch(e_em)
#import patches 
qbox = patches.Rectangle((4.5-3.2,-4),2.5,8,linewidth=1,edgecolor='k',facecolor='xkcd:blue grey',alpha=0.2) 
ax1 = plt.gca()
#ax1.add_patch(qbox)
#plt.text(4.5-2.4,-0.1,'???????',color='k',fontsize=20)
#min_decay = np.str(decay_const.min()*100)
#max_decay = np.str(decay_const.max()*100)
#rate = np.str(decay_const*100)
plt.legend([lgd1,lgd2,lgd3,lgd4,lgd5],['Sluggish Archean Plate textonics','Previous estimates','This study','Early Archean emergence','Late Archean emergence'],\
           loc='upper left',bbox_to_anchor=(0.1, 1),fontsize=12)#,'Late continent emergence'
plt.plot(lear_age_Ga,lear_d18O,color='xkcd:grey blue')

plt.ylim([-2.5,7])
plt.xlim([0,4.5])

inset1 = fig3.add_axes([.6, .52, .32, .25])

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}

matplotlib.rc('font', **font)

lear_ticks = ['60','40','20','0']
inset1.plot(lear_age,lear_d18O,color='xkcd:grey blue')
plt.errorbar(np.multiply(1e3,inversion_ages[2:]),inversion_d18O[2:],yerr=inversion_err[2:],fmt='o',marker='o', mfc=this_study_color,\
         mec=this_study_color,ecolor=this_study_color,capsize=10, elinewidth=2)
plt.xlabel('Age (Ma)',fontsize=10);plt.ylabel('Seawater $\delta^{18}$O (‰ VSMOW)',fontsize=10)
ins_ax = plt.gca()
ins_ax.set_yticks([-1,-0.5,0])
plt.legend(['Benthic foraminifera record'],fontsize=10)


ax3=plt.gca()
ax3.invert_xaxis()
#ax2.set_aspect(0.2)
fig3.set_size_inches(12, 8.5)    
fig3.savefig('ocean_d18O.eps',metadata='eps')
 



#%%
#### Supplemental figures
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)

num_rows = 4; num_cols = 2
emergence_fig= plt.figure(); plt.clf()

aa=np.exp(-np.multiply(time,decay_const_low*k_sum[-1]))
bb=1-aa 

k_list = [k_weath,k_growth,k_hiT,k_loT,k_W_recycling]
frac_flux = np.zeros(len(k_list))
for ii in range(0,len(k_list)): 
    frac_flux[ii]=k_list[ii]/k_sum[-1]
    
kdecay = np.zeros([len(k_list),time.size])


for ientry in range(0,len(k_list)):
    for iientry in range(0,len(time)):
        if time[iientry]<1.5:
            kdecay[ientry,iientry] = decay_const_low*k_list[ientry]
        else:
            kdecay[ientry,iientry] = k_list[ientry]


ocean_scenarios = [dW_decay_low,dW_early,dW_middle,dW_late]

forcings = [kdecay,
            [k_weathering_early,k_growth_early,k_hiT_change,k_loT_change,k_water_change],
            [k_weathering_change,k_growth_change,k_hiT_change,k_loT_change,k_water_change],
            [k_weathering_late,k_growth_late,k_hiT_change,k_loT_change,k_water_change]]

idx = -1
symb_list = ['--',':','-.','-','-']
marker_list = ['x','','','.','']
rate = np.str(decay_const_low*100)
title_list = ['Sluggish Archean plate tectonics','Early Archean emergence','Late Archean emergence',
              'Proterozoic continental emergence']
legend_loc = [0.1,2.2]
time_labels = ['4.5','3.5','2.5','1.5','0.5']
time_ticks = np.linspace(0,4,5)
for irow in range(0,len(ocean_scenarios)):
    idx=idx+num_cols
    a1=plt.subplot(num_rows,num_cols,idx)
    a1.errorbar(np.subtract(4.5,previous_ages),previous_d18O,yerr=previous_err,fmt='o',marker='o', mfc='xkcd:silver',\
         mec='xkcd:silver',ecolor='xkcd:silver',capsize=10, elinewidth=20)

    a1.errorbar(np.subtract(4.5,inversion_ages[1:]),inversion_d18O[1:],yerr=inversion_err[1:],fmt='o',marker='o', mfc=this_study_color,\
         mec=this_study_color,ecolor=this_study_color,capsize=10, elinewidth=2)

    a1.errorbar(np.subtract(4.5,inversion_ages[0]),inversion_d18O[0],yerr=inversion_err[0],fmt='o',marker='o', mfc=this_study_color,\
         mec=this_study_color,ecolor=this_study_color,capsize=10, elinewidth=20)
    
    a1.plot(time_new,ocean_scenarios[irow],'k')
    a1.set_title(title_list[irow],x=1.1)
    a1.set_xticks(time_ticks)
    a2=plt.subplot(num_rows,num_cols,idx+1)
    a2.set_xticks(time_ticks)
    a2.yaxis.set_label_position("right")
    a2.yaxis.set_ticks_position("right")
    for ientry in range(0,np.size(forcings,axis=1)):
        a2.plot(time,forcings[irow][ientry],'k',marker=marker_list[ientry],linestyle=symb_list[ientry])
    
    if irow==0:
        a1.legend(['Seawater $\delta^{18}$O','Previous constraints','This study'],\
           loc='upper left',bbox_to_anchor=legend_loc,fontsize=10)
        lgd=a2.legend(['Cont. weathering','Continental Recycling','High T alteration','Low T alteration','Water Recycling'],\
           loc='upper left',bbox_to_anchor=legend_loc,fontsize=10)
    if irow==3:
        a1.set_xlabel('Age (Ga)')
        a1.set_xticklabels(time_labels)
        a2.set_xlabel('Age (Ga)')
        a2.set_xticklabels(time_labels)
    else: 
        a1.set_xticklabels('')
        a2.set_xticklabels('')
    if irow==2:
        a1.set_ylabel('Seawater $\delta^{18}$O (‰ VSMOW)')
        a2.set_ylabel('Rate (Gyr $^{-1}$)')
    
emergence_fig.subplots_adjust(wspace=0.08,hspace=0.5)    
emergence_fig.set_size_inches(12, 9.5)    
emergence_fig.savefig('emergence_scenarios.eps',bbox_extra_artists=(lgd,),bbox_inches='tight',metadata='eps')

#%%
fig4, (ax1, ax2) = plt.subplots(1, 2, figsize=[5.5, 2.8]) 
binwidth = .1
ax1.hist(wing_W_R,bins=np.arange(min(wing_W_R), max(wing_W_R) + binwidth, binwidth), weights=wing_weights,color=wing_color,ec='k',alpha=0.5)
ax1.hist(hoku50_W_R,bins=np.arange(min(hoku50_W_R), max(hoku50_W_R) + binwidth, binwidth), weights=hoku50_weights,color=hoku50_color,ec='k',alpha=0.5)
ax1.hist(skaer_W_R,bins=np.arange(min(skaer_W_R), max(skaer_W_R) + binwidth, binwidth), weights=skaer_weights,color=skaer_color,ec='k',alpha=0.5)
ax1.set_title('Forward models')
ax1.set_xlabel('F/R ratio'); ax1.set_ylabel('Fraction of runs')
ax1.legend(['Equilibrium','Kinetic','Skaergaard'],loc='upper center')


binwidth = .1
ax2.hist(ohmoto_W_R,bins=np.arange(min(ohmoto_W_R), max(ohmoto_W_R) + binwidth, binwidth), weights=ohmoto_weights,color=ohmoto_color,ec='k',alpha=alphabar)
ax2.hist(EPR_W_R,bins=np.arange(min(EPR_W_R), max(EPR_W_R) + binwidth, binwidth), weights=EPR_weights,color=EPR_color,ec='k',alpha=alphabar)
ax2.hist(solea_W_R,bins=np.arange(min(solea_W_R), max(solea_W_R) + binwidth, binwidth), weights=solea_weights,color=solea_color,ec='k',alpha=alphabar)
ax2.hist(Whole_W_R,bins=np.arange(min(Whole_W_R), max(Whole_W_R) + binwidth, binwidth), weights=pano_weights,color=pano_color,ec='k',alpha=alphabar)
ax2.set_title('Crustal Sections')
#ax2 = plt.gca()
ax2.set_yticklabels('')
ax2.legend(['Fukazawa','EPR','Solea','Panorama'],loc='upper right')#,bbox_to_anchor=(-0.1, 1)
ax2.set_xlabel('F/R ratio'); 

fig4.set_size_inches(12, 8.5)    
fig4.savefig('W_R_results.eps',metadata='eps')

#%%
fig_tempiso = plt.figure(); plt.clf()
iso_grid_list = [wing_iso_grid,hoku50_iso_grid,skaer_iso_grid,EPR_iso_grid,ohmoto_iso_grid,solea_iso_grid,
                 Whole_iso_grid]
#iso_grid_list = [EPR_iso_grid,ohmoto_iso_grid,solea_iso_grid] #modern settings only
temp_grid_list = [temp_grid_wing,temp_grid_hoku50,temp_grid_skaer,temp_grid_EPR,temp_grid_ohmoto,temp_grid_solea,
                  temp_grid_Whole]

#temp_grid_list = [temp_grid_EPR,temp_grid_ohmoto,temp_grid_solea]

x_list = np.divide([x_wing,x_hoku50,x_skaer/1000,x_EPR,x_ohmoto,x_solea,X_Whole],1000)
#x_list = np.divide([x_EPR,x_ohmoto,x_solea],1000)

z_list = np.divide([z_wing,z_hoku50,z_skaer,z_EPR,z_ohmoto,z_solea,Z_Whole],1000)
#z_list = np.divide([z_EPR/10,z_ohmoto,z_solea],1000)

num_yticks = 3
z_ticks = np.zeros([len(z_list),num_yticks])
for i in range(len(z_list)):
    z_ticks[i] = np.around(np.linspace(z_list[i].min(),z_list[i].max(),num_yticks),decimals=1)
    
idx=-2
title_list = ['Equilibrium','Kinetic','Skaer','EPR','Fukazawa','Solea','Pano']
#title_list = ['EPR','Fukazawa','Solea']

from matplotlib.ticker import LinearLocator
for i in range(len(iso_grid_list)):
    idx=idx+3
#    print(idx)
    a1 = plt.subplot(7,3,idx)
    a1con = a1.contourf(iso_grid_list[i],cmap=isogrid_cmap, extent=[0-0.5, temp_grid_list[i].max(),0-0.5, temp_grid_list[i].max()])
    
    numticks = len(a1.get_xticks())
    xlims = a1.get_xlim(); ylims=a1.get_ylim()
    a1.set_xticklabels(np.around(np.linspace(x_list[i].min(),x_list[i].max(),numticks+1)))
    a1.get_yaxis().set_major_locator(LinearLocator(numticks=num_yticks))
    a1.set_yticklabels(z_ticks[i])
    
    
    tickrange = np.round(np.linspace(iso_grid_list[i].min(),iso_grid_list[i].max(),num=3),decimals=1)
    cbar = plt.colorbar(a1con)
    cbar.set_ticks(tickrange)

       
    a2 = plt.subplot(7,3,idx+1)
#    a2con = a2.contourf(xi,yi,zii,levels,cmap=tempgrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
    a2con = a2.contourf(temp_grid_list[i],cmap=tempgrid_cmap, extent=[0-0.5, temp_grid_list[i].max(),0-0.5, temp_grid_list[i].max()])
    a2.set_xticklabels(np.around(np.linspace(x_list[i].min(),x_list[i].max(),numticks+1)))
    cbar2 = plt.colorbar(a2con)
    tempticks = np.around(np.linspace(temp_grid_list[i].min(),temp_grid_list[i].max(),num=3),decimals=0)
    cbar2.set_ticks(tempticks)
    plt.ylabel('');a2.set_yticklabels('')
    
    
    a3 = plt.subplot(7,3,idx+2)
    a3.scatter(x_list[i],z_list[i],marker='.')
    a3.set_yticklabels('');a3.set_xticklabels('')
#    a3.text(0.7*xlims[1],ylims[0]+0.2*ylims[1],title_list[i],fontsize=10,color='w',backgroundcolor='k')
    a3.set_title(title_list[i],fontsize=11,color='k') #,backgroundcolor='k' ,y=1.02
#    a3.set_
    a3.set_aspect('equal')
    
    if idx==1:
        a1.set_title(r'$\delta^{18}$O (‰ VSMOW)')
        a2.set_title('Temp. ($^\circ$C)')
        
    if idx==19: #19
        a1.set_xlabel('Horizontal (km)')
        a2.set_xlabel('Horizontal (km)')
        a3.set_xlabel('True Aspect Ratio')
        a1.contourf(xi,yi,zi,cmap=isogrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
        a2.contourf(xi,yi,zii,levels,cmap=tempgrid_cmap,extent=[0-0.5, X_Whole.max(),0-0.5, Z_Whole.max()])
#        a3.set_xlabel('Horizontal (km)')
    
    if idx==10: #10
        a1.set_ylabel('Vertical (km)')
#    del(a1,a2,a1con,a2con)

plt.subplots_adjust(wspace=0.15,hspace=.7)
fig_tempiso.set_size_inches(14,8.5)
fig_tempiso.savefig('temp_iso_all.eps',metadata='eps')

#%%
