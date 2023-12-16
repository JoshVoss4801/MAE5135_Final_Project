import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from FlowPy_LSM_Main import *

#===============[ Gravity ]===============#
g = 9.81 # Gravitational constant
g_xDir = 0 # Gravity in x-Direction
g_yDir = 0 # Gravity in y-Direction
grav = Gravity(g,g_xDir,g_yDir)
#########################[ Grid Parameters ]#########################
LC = 0.02 # Characteristic length = Initial Bubble Radius
LX = 0.06 # x length
LY = 0.06 # y height
NX = 201 # Number of x steps !KEEP ODD!
NY = 201 # Number of y steps !KEEP ODD!
cavity = Space(NX,NY,LX,LY,LC)
#########################[ Define Initial Constants ]#########################
Re = 100 # Reynolds Number
We = 30 # Weber Number
#===============[ Fluid Properties ]===============# 
rho_r = 100 # Ratio of rho1/rho1 where rho1 is the density of the fluid outside the interface
mu_r = 100 # Ratio of mu1/mu1 where mu1 is the veiscosity of the fluid outside the interface
sigma = 0.072 #dont think its needed?
#===============[Time Parameters]===============#
time = 0.1
dt = 0.01

#mu1 = 18.13e-6
mu1 = 1e-3
rho1 = 1e3
GetInitialConstants(cavity,grav,Re,We,rho_r,mu_r,rho1,mu1,sigma,dt)

#########################[ Define Initial Level Set Function ]#########################
r = 0.02
x0 = 0.03
y0 = 0.03
a = 1
b = 1
GetInitialPhiBubbleLoop(cavity,x0,y0,a,b,r) # Define Initial Phi
GetH(cavity)
GetInitialPressure(cavity)
cavity.SetCentre()
PlotPressureContour(cavity)

#########################[ Begin Flow Simulation ]#########################
print("######## Beginning FlowPy Simulation ########")
print("#############################################")
print("# Simulation time: {0:.2f}".format(time))
#print("# Mesh: {0} x {1}".format(NY,NX))
#print("# Re/u: {0:.2f}\tRe/v:{1:.2f}".format(cavity.Re,cavity.Re))

t = 0
i = 0
# Time Loop
while(t<time):
    # Print time left
    #sys.stdout.write("\rSimulation time left: {0:.2f}".format(time-t))
    #sys.stdout.flush()
    print(f"Time Step: {i}")

    SetBoundaryCondition(cavity)

    GetDelta(cavity)
    GetCurvature(cavity)
    GetH(cavity)

    GetDensity(cavity)
    GetViscosity(cavity)

    cavity.SetCentre()
    PlotQuad(cavity,dt)
    
    GetStarredVelocitiesNEW(cavity)
    SolvePressurePoisson(cavity)
    SolveMomentumEquation(cavity)

    GetNewLSF(cavity)
    
    t+=dt
    i+=1


