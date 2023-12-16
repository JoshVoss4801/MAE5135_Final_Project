import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import os

class Gravity:
    def __init__(self,grav_const,g_xDir,g_yDir):
        self.GetGravity(grav_const,g_xDir,g_yDir)
    
    def GetGravity(self,grav_const,g_xDir,g_yDir):
        self.g = float(grav_const)
        self.gx = float(g_xDir)
        self.gy = float(g_yDir)
class Space:
    def __init__(self,colpts,rowpts,length,breadth,char_length):
        self.CreateMesh(rowpts,colpts)
        self.SetDeltas(breadth,length,char_length)

    def CreateMesh(self,rowpts,colpts):
        #Domain gridpoints
        self.rowpts = rowpts
        self.colpts = colpts

        #Main matrices
        self.u = np.zeros((self.rowpts+2,self.colpts+2)) # Velocity x-Direction
        self.v = np.zeros((self.rowpts+2,self.colpts+2)) # Velocity y-Direction
        self.u_star = np.zeros((self.rowpts+2,self.colpts+2)) # u_star
        self.v_star = np.zeros((self.rowpts+2,self.colpts+2)) # v_star
        self.p = np.zeros((self.rowpts+2,self.colpts+2)) # Pressure - x
        self.phi = np.zeros((self.rowpts+2,self.colpts+2)) # Level set function value
        self.rho = np.zeros((self.rowpts+2,self.colpts+2)) # Density
        self.mu = np.zeros((self.rowpts+2,self.colpts+2)) # Viscosity
        self.H = np.zeros((self.rowpts+2,self.colpts+2)) # H for interface
        self.k = np.zeros((self.rowpts+2,self.colpts+2)) # Kurvature
        self.delta = np.zeros((self.rowpts+2,self.colpts+2)) #delta fuction for interface

        #Centered Matrixs used for plotting and output
        self.u_c = np.zeros((self.rowpts,self.colpts))
        self.v_c = np.zeros((self.rowpts,self.colpts))
        self.u_star_c = np.zeros((self.rowpts,self.colpts))
        self.v_star_c = np.zeros((self.rowpts,self.colpts))
        self.p_c = np.zeros((self.rowpts,self.colpts))
        self.phi_c = np.zeros((self.rowpts,self.colpts))
        self.rho_c = np.zeros((self.rowpts,self.colpts))
        self.mu_c = np.zeros((self.rowpts,self.colpts))
        self.H_c = np.zeros((self.rowpts,self.colpts))
        self.k_c = np.zeros((self.rowpts,self.colpts))
        self.delta_c = np.zeros((self.rowpts,self.colpts))
    def SetDeltas(self,breadth,length,char_length):
        lc = float(char_length)
        lx = float(length)
        ly = float(breadth)

        dx = lx/(self.colpts-1)
        dy = ly/(self.rowpts-1)
        epsilon = dx+dy


        self.lc = float(char_length)
        self.lx = float(length)
        self.ly = float(breadth)
        self.dx = float(dx)
        self.dy = float(dy)
        self.epsilon = float(epsilon)
    def SetCentre(self):
        self.u_c=self.u[1:-1,1:-1]
        self.v_c=self.v[1:-1,1:-1]
        self.u_star_c=self.u_star[1:-1,1:-1]
        self.v_star_c=self.v_star[1:-1,1:-1]
        self.p_c=self.p[1:-1,1:-1]
        self.phi_c=self.phi[1:-1,1:-1]
        self.rho_c=self.rho[1:-1,1:-1]
        self.mu_c=self.rho[1:-1,1:-1]
        self.H_c=self.H[1:-1,1:-1]
        self.k_c=self.H[1:-1,1:-1]
        self.delta_c=self.H[1:-1,1:-1]


def GetInitialConstants(space,gravity,Re,We,rhor,mur,rho1,mu1,sigma,dt):
    lc = float(space.lc)
    re = float(Re)
    we = float(We)
    muo = float(mu1)
    sig = float(sigma)
    g = gravity.g

    uc = (sig/muo)*(We/Re)
    Fr = uc/np.sqrt(g*lc)

    space.dt = float(dt)
    space.rho_r = rhor
    space.mu_r = mur
    space.rho1 = float(rho1)
    space.mu1 = float(mu1)
    space.sigma = sigma
    space.Re = re
    space.We = we
    space.Fr = Fr
    space.gx = float(gravity.gx)
    space.gy = float(gravity.gy)
    space.uc = uc

    
    
def GetInitialPhiBubbleLoop(space,x0,y0,a,b,r):
    phi = space.phi.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dy = float(space.dy)
    dx = float(space.dx)
    rad = float(r)
    X0 = float(x0)
    Y0 = float(y0)
    A = float(a)
    B = float(b)
    

    y = np.zeros(rows+2)
    x = np.zeros(cols+2)

    #set matrix space for initiap phi
    intphi = phi.copy()

    for j in range(0,rows+2):
        y[j] = (2*j-1)*dy/2
        dy_sq = abs((y[j]-Y0)**2/B)
        for i in range(0,cols+2):
            x[i] = (2*i-1)*dx/2
            dx_sq = abs((x[i]-X0)**2/A)
            #phi[j,i] = np.sqrt((x[i]-x0)**2/a + (y[j]-y0)**2/b) - r
            #phi[j,i] = ((x[i]-x0)*(x[i]-x0)/a + (y[j]-y0)*(y[j]-y0)/b)**(1/2) - r
            intphi[j,i] = np.sqrt(dx_sq+dy_sq) - rad

    space.phi = intphi.copy()
def GetInitialPhiBubble(space,x0,y0,a,b,r):
    phi = space.phi.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dy = float(space.dy)
    dx = float(space.dx)
    rad = float(r)
    X0 = float(x0)
    Y0 = float(y0)
    A = float(a)
    B = float(b)
def SetTimeStep_Test(Dt,space):
    dt = Dt
    space.dt = dt

def SetBoundaryCondition(space):
    u = space.u.astype(float,copy=False)
    v = space.v.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)

    space.u[:,0] = u[:,1]
    space.u[:,cols+1] = u[:,cols]
    space.u[0,:] = u[1,:]
    space.u[rows+1,:] = u[rows,:]

    space.v[:,0] = v[:,1]
    space.v[:,cols+1] = v[:,cols]
    space.v[0,:] = v[1,:]
    space.v[rows+1,:] = v[rows,:]
def GetInitialPressure(space):
    rows = int(space.rowpts)
    cols = int(space.colpts)
    rho1 = float(space.rho1)
    mu1 = float(space.mu1)
    Re = float(space.Re)
    We = float(space.We)

    H = space.H.astype(float,copy=False)
    p = space.p.astype(float,copy=False)


    for j in range(0,rows+2):
        for i in range(0,cols+2):
            if H[j,i] == 0:
                p[j,i] = 2*mu1*Re/(rho1*We)
            elif H[j,i] == 1:
                p[j,i] = 0
            else:
                p[j,i] = (1-H[j,i])*(2*mu1*Re/(rho1*We))

    space.p = p.copy()
def GetH(space):
    phi = space.phi.astype(float,copy=False)
    H = space.H.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    e = float(space.epsilon)

    for j in range(0,rows+2):
        for i in range(0,cols+2):
            if phi[j,i] < -e:
                H[j,i] = 0
            elif -e <= phi[j,i] and phi[j,i] <= e:
                H[j,i] = (phi[j,i]+e)/2/e + np.sin(np.pi*phi[j,i]/e)/2/np.pi
            else:
                H[j,i] = 1

    space.H = H.copy()
def GetDelta(space):
    phi = space.phi.astype(float)
    delta = space.delta.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    e = float(space.epsilon)

    for j in range(0,rows+2):
        for i in range(0,cols+2):
            if -e > phi[j,i]:
                delta[j,i] = 0
            elif abs(phi[j,i]) <= e:
                delta[j,i] = 1/(2*e) + np.cos(np.pi*phi[j,i]/e)/(2*e)
            else:
                delta[j,i] = 0

    space.delta = delta.copy()
def GetCurvature(space):
    k = space.k.astype(float,copy=False)
    phi = space.phi.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dx = float(space.dx)
    dy = float(space.dy)

    # Calculate differential Terms
    phi1_x = (phi[1:rows+1,2:]-phi[1:rows+1,0:cols])/(2*dx) 
    phi1_y = (phi[2:,1:cols+1]-phi[0:rows,1:cols+1])/(2*dy)

    phi2_x = (phi[1:rows+1,2:]-2*phi[1:rows+1,1:cols+1]+phi[1:rows+1,0:cols])/dx**2
    phi2_y = (phi[2:,1:cols+1]-2*phi[1:rows+1,1:cols+1]+phi[0:rows,1:cols+1])/dy**2
    phi2_xy = ( (phi[2:,2:]-phi[0:rows,2:])/(2*dy) - (phi[2:,0:cols]-phi[0:rows,0:cols])/(2*dy) )/(2*dx)

    k[1:rows+1,1:cols+1] = (phi2_x + phi2_y)/np.sqrt(phi1_x**2 + phi1_y**2) \
                            - ( phi1_x*(phi1_x*phi2_x + phi1_y*phi2_xy) + phi1_y*(phi1_x*phi2_xy + phi1_y*phi2_y) )/((phi1_x**2 + phi1_y**2)**(3/2))

    space.k = k.copy()

def GetDensity(space):
    H = space.H.astype(float,copy=False)
    #rho1 = float(space.rho1)
    rho_r = float(space.rho_r)
    #rho = space.rho.astype(float,copy=False)

    rho = 1 + (rho_r-1)*H
    
    space.rho = rho.copy()
def GetViscosity(space):
    H = space.H.astype(float,copy=False)
    #mu1 = float(space.mu1)
    mu_r = float(space.mu_r)
    #rho = space.rho.astype(float,copy=False)

    mu = 1 + (mu_r-1)*H
    
    space.mu = mu.copy()

def GetStarredVelocities(space):
    u = space.u.astype(float,copy=False)
    v = space.v.astype(float,copy=False)
    phi = space.phi.astype(float,copy=False)
    rho = space.rho.astype(float,copy=False)
    mu = space.mu.astype(float,copy=False)
    k = space.k.astype(float,copy=False)
    delta = space.delta.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    Re = float(space.Re)
    Fr = float(space.Fr)
    We = float(space.We)
    lc = float(space.lc)
    gx = float(space.gx)
    gy = float(space.gy)
    sigma = float(space.sigma)

    # Calculate differential terms
    u1_x = (u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx) # du/dx, (Ui+1 - Ui-1)/2dx
    u1_y = (u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy) # du/dy, (Uj+1 - Uj-1)/2dy
    v1_x = (v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx) # dv/dx, (Vi+1 - Vi-1)/2dx
    v1_y = (v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy) # dv/dy, (Vj+1 - Vj-1)/2dy
    phi1_x = (phi[1:rows+1,2:]-phi[1:rows+1,0:cols])/(2*dx) # dphi/dx, (PHIi+1 - PHIi-1)/2dx
    phi1_y = (phi[2:,1:cols+1]-phi[0:rows,1:cols+1])/(2*dy) # dphi/dy, (PHIi+1 - PHIi-1)/2dy

    u2_x = (u[1:rows+1,2:]-2*u[1:rows+1,1:cols+1]+u[1:rows+1,0:cols])/dx**2 # ddu/dxdx, (Ui+1 -  2U + Ui-1)/dx^2
    u2_y = (u[2:,1:cols+1]-2*u[1:rows+1,1:cols+1]+u[0:rows,1:cols+1])/dy**2 # ddu/dydy, (Uj+1 -  2U + Uj-1)/dy^2
    v2_x = (v[1:rows+1,2:]-2*v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols])/dx**2 # ddv/dxdx, (Vi+1 -  2V + Vi-1)/dx^2
    v2_y = (v[2:,1:cols+1]-2*v[1:rows+1,1:cols+1]+v[0:rows,1:cols+1])/dy**2 # ddv/dydy, (Vj+1 -  2V + Vj-1)/dy^2

    conv_x = u[1:rows+1,1:cols+1]*u1_x + v[1:rows+1,1:cols+1]*u1_y
    diff_x = (mu[1:rows+1,1:cols+1]/(rho[1:rows+1,1:cols+1]*Re))*(u2_x+u2_y)
    grav_x = gx/Fr
    Bx = sigma*k[1:rows+1,1:cols+1]*delta[1:rows+1,1:cols+1]*phi1_x*lc**2
    surf_x = Bx/(We*rho[1:rows+1,1:cols+1])

    conv_y = u[1:rows+1,1:cols+1]*v1_x + v[1:rows+1,1:cols+1]*v1_y
    diff_y = (mu[1:rows+1,1:cols+1]/(rho[1:rows+1,1:cols+1]*Re))*(v2_x+v2_y)
    grav_y = gy/Fr
    By = sigma*k[1:rows+1,1:cols+1]*delta[1:rows+1,1:cols+1]*phi1_y*lc**2
    surf_y = By/(We*rho[1:rows+1,1:cols+1])

    print(f"dt*conv_x: {np.amax(dt*conv_x)}")
    print(f"dt*diff_x: {np.amax(dt*diff_x)}")
    print(f"dt*grav_x: {np.amax(dt*grav_x)}")
    print(f"dt*surf_x: {np.amax(dt*surf_x)}")
    print(f"dt*conv_y: {np.amax(dt*conv_y)}")
    print(f"dt*diff_y: {np.amax(dt*diff_y)}")
    print(f"dt*grav_y: {np.amax(dt*grav_y)}")
    print(f"dt*surf_y: {np.amax(dt*surf_y)}")

    u_star = u.copy()
    v_star = v.copy()

    u_star[1:rows+1,1:cols+1] = u[1:rows+1,1:cols+1] + dt*(-conv_x+diff_x+grav_x+surf_x)
    v_star[1:rows+1,1:cols+1] = v[1:rows+1,1:cols+1] + dt*(-conv_y+diff_y+grav_y+surf_y)

    space.u_star = u_star.copy()
    space.v_star = v_star.copy()
def SolvePressurePoisson(space):
    u_star = space.u_star.astype(float,copy=False)
    v_star = space.v_star.astype(float,copy=False)
    p = space.p.astype(float,copy=False)
    rho = space.rho.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dt = float(space.dt)
    dx = float(space.dx)
    dy = float(space.dy)
    lx = float(space.lx)
    ly = float(space.ly)

    #calculate cell density values
    rhoE = (rho[1:rows+1,2:]+rho[1:rows+1,1:cols+1])/2
    rhoW = (rho[1:rows+1,0:cols]+rho[1:rows+1,1:cols+1])/2
    rhoN = (rho[2:,1:cols+1]+rho[1:rows+1,1:cols+1])/2
    rhoS = (rho[0:rows,1:cols+1]+rho[1:rows+1,1:cols+1])/2
    
    #Denominator
    den = ((1/rhoE)+(1/rhoW))/dx**2 + (1/rhoN)+(1/rhoS)/dy**2

    # Calculate differential terms
    ustar1_x = (u_star[1:rows+1,2:]-u_star[1:rows+1,0:cols])/(2*dx)
    vstar1_y = (v_star[2:,1:cols+1]-v_star[0:rows,1:cols+1])/(2*dy)

    error = 1
    tol = 1e-3
    ########################################
    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    ########################################
    i = 0
    while (error>tol):
        i+= 1

        # Save current pressure as p_temp
        p_temp = p.astype(float,copy=True)

        # P_TEMP Terms
        p_temp_terms = ((p_temp[1:rows+1,2:]/rhoE)+(p_temp[1:rows+1,0:cols]/rhoW))/dx**2 + ((p_temp[2:,1:cols+1]/rhoN)+(p_temp[0:rows,1:cols+1]/rhoS))/dy**2
        # VEL Terms
        vel_terms = (ustar1_x+vstar1_y)

        p[1:rows+1,1:cols+1] = (vel_terms+p_temp_terms)/den

        
        p_temp_c = p_temp[1:-1,1:-1]
        p_c = p[1:-1,1:-1]
        plt.subplot(1, 2, 1)
        plt.colorbar(plt.contourf(X, Y, p_temp_c))
        plt.subplot(1, 2, 2)
        plt.colorbar(plt.contourf(X, Y, p_c))
        plt.tight_layout()
        #plt.show()
        # Apply pressure boundary conditions
        #SetPBoundary(space,left,right,top,bottom)
        p[:,0] = p[:,1] # Left
        p[:,cols+1] = p[:,cols]
        p[0,:] = p[1,:]
        p[rows+1,:] = p[rows,:]

        # Find maximum error between old and new pressure matrices
        error=np.amax(abs(p-p_temp))
        
        #Escape condition in case solution does not converge after 500 iterations
        if(i>500):
            tol*=10
            print(f"error: {error}")
    
    space.p = p.copy()
def SolveMomentumEquation(space):
    u = space.u.astype(float,copy=False)
    v = space.v.astype(float,copy=False)
    rho = space.rho.astype(float,copy=False)
    p = space.p.astype(float,copy=False)
    u_star = space.u_star.astype(float,copy=False)
    v_star = space.v_star.astype(float,copy=False)
    
    rows = int(space.rowpts)
    cols = int(space.colpts)
    dt = float(space.dt)
    dx = float(space.dx)
    dy = float(space.dy)

    # Calculate pressure diffentials
    p1_x = (p[1:rows+1,2:]-p[1:rows+1,0:cols])/(2*dx)
    p1_y = (p[2:,1:cols+1]-p[0:rows,1:cols+1])/(2*dy)
    print(f"dt*px: {np.amax(dt*p1_x/rho[1:rows+1,1:cols+1])}")
    print(f"dt*py: {np.amax(dt*p1_y/rho[1:rows+1,1:cols+1])}")

    #Calculate u and v at next timestep
    u[1:rows+1,1:cols+1] = u_star[1:rows+1,1:cols+1] - dt*p1_x/rho[1:rows+1,1:cols+1]
    v[1:rows+1,1:cols+1] = v_star[1:rows+1,1:cols+1] - dt*p1_y/rho[1:rows+1,1:cols+1]

    space.u = u.copy()
    space.v = v.copy()

def GetNewLSF(space):
    u = space.u.astype(float,copy=False)
    v = space.v.astype(float,copy=False)
    phi = space.phi.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    
    phi1_x = (phi[1:rows+1,2:]-phi[1:rows+1,0:cols])/(2*dx)
    phi1_y = (phi[2:,1:cols+1]-phi[0:rows,1:cols+1])/(2*dy)
    phi[1:rows+1,1:cols+1] = phi[1:rows+1,1:cols+1] - dt*(u[1:rows+1,1:cols+1]*phi1_x + v[1:rows+1,1:cols+1]*phi1_y)

    space.phi = phi.copy()

def PlotPressureContour(space):
    p_p = space.p_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,p_p)
    fig.colorbar(cont)
    fig.tight_layout()
    fig.show()
def PlotDensityContour(space):
    rho_p = space.rho_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,rho_p)
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()
def PlotHContour(space):
    H_p = space.H_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,H_p)
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()
def PlotPhiContour(space):
    phi_p = space.phi_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,phi_p)
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()
def PlotKContour(space):
    k_p = space.k_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,k_p)
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()
def PlotDeltaContour(space):
    delta_p = space.delta_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    #index_cut_x=int(cols/10)
    #index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,delta_p)
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()
def PlotVel(space):
    u_p = space.u_c.astype(float,copy=False)
    v_p = space.v_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)


    #magVel_p = np.sqrt(u_p**2 + v_p**2)

    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    #Determine indexing for stream plot (10 points only)
    index_cut_x=int(cols/10)
    index_cut_y=int(rows/10)
    #Create blank figure
    fig=plt.figure(figsize=(10,10))
    ax=plt.axes(xlim=(0,lx),ylim=(0,ly))
    #Create initial contour and stream plot as well as color bar
    ax.set_xlim([0,lx])
    ax.set_ylim([0,ly])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    cont=ax.contourf(X,Y,u_p)
    stream = ax.streamplot(X[::index_cut_y,::index_cut_x],Y[::index_cut_y,::index_cut_x],u_p[::index_cut_y,::index_cut_x],v_p[::index_cut_y,::index_cut_x],color="k")
    fig.colorbar(cont)
    fig.tight_layout()
    plt.show()

def PlotQuad(space,timestep):
    u_p = space.u_c.astype(float,copy=False)
    v_p = space.v_c.astype(float,copy=False)
    u_star_p = space.u_star_c.astype(float,copy=False)
    v_star_p = space.v_star_c.astype(float,copy=False)
    phi_p = space.phi_c.astype(float,copy=False)
    rho_p = space.rho_c.astype(float,copy=False)
    p_p = space.p_c.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    lx = float(space.lx)
    ly = float(space.ly)

    x=np.linspace(0,lx,cols)
    y=np.linspace(0,ly,rows)
    [X,Y]=np.meshgrid(x,y)
    
    index_cut_x=int(cols/10)
    index_cut_y=int(rows/10)

    # Create a larger figure
    fig = plt.figure(figsize=(10, 10))

    # Create a 2x2 subplot grid
    plt.subplot(2, 2, 1)  # Subplot in the top-left position
    #plt.axes(xlim=(0,lx),ylim=(0,ly))
    plt.colorbar(plt.contourf(X, Y, u_p))
    plt.streamplot(X[::index_cut_y,::index_cut_x],Y[::index_cut_y,::index_cut_x],u_p[::index_cut_y,::index_cut_x],v_p[::index_cut_y,::index_cut_x],color="k")
    plt.title('u')

    plt.subplot(2, 2, 2)  # Subplot in the top-right position
    #plt.axes(xlim=(0,lx),ylim=(0,ly))
    plt.colorbar(plt.contourf(X, Y, u_star_p))
    plt.streamplot(X[::index_cut_y,::index_cut_x],Y[::index_cut_y,::index_cut_x],u_star_p[::index_cut_y,::index_cut_x],v_star_p[::index_cut_y,::index_cut_x],color="k")
    plt.title('u_star')

    plt.subplot(2, 2, 3)  # Subplot in the bottom-left position
    #plt.axes(xlim=(0,lx),ylim=(0,ly))
    plt.colorbar(plt.contourf(X, Y, rho_p))
    plt.title('rho')

    plt.subplot(2, 2, 4)  # Subplot in the bottom-right position
    #plt.axes(xlim=(0,lx),ylim=(0,ly))
    plt.colorbar(plt.contourf(X, Y, p_p))
    plt.title('p')

    # Add a title to the entire figure including the timestep
    plt.suptitle(f'Timestep: {timestep}', fontsize=16)

    # Adjust layout for better spacing
    plt.tight_layout()

    # Show the figure
    plt.show()

def GetStarredVelocitiesNEW(space):

    u = space.u.astype(float,copy=False)
    v = space.v.astype(float,copy=False)
    phi = space.phi.astype(float,copy=False)
    rho = space.rho.astype(float,copy=False)
    mu = space.mu.astype(float,copy=False)
    k = space.k.astype(float,copy=False)
    delta = space.delta.astype(float,copy=False)

    rows = int(space.rowpts)
    cols = int(space.colpts)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    lc = float(space.lc)
    Re = float(space.Re)
    We = float(space.We)
    Fr = float(space.Fr)
    sigma = float(space.sigma)
    gx = float(space.gx)
    gy = float(space.gy)

    # Central differentials
    dphidx = (phi[1:rows+1,2:]-phi[1:rows+1,0:cols])/(2*dx) # dphi/dx
    dphidy = (phi[2:,1:cols+1]-phi[0:rows,1:cols+1])/(2*dy) # dphi/dy
    dudx = (u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx) # du/dx
    dudy = (u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy) # du/dy
    dvdx = (v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx) # dv/dx
    dvdy = (v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy) # dv/dy

    # Cell Wall Differentials
    dudxE = (u[1:rows+1,2:]-u[1:rows+1,1:cols+1])/dx # du/dx on East cell wall
    dudxW = (u[1:rows+1,1:cols+1]-u[1:rows+1,0:cols])/dx # du/dx on West cell wall
    dudyN = (u[2:,1:cols+1]-u[1:rows+1,1:cols+1])/dy # du/dy on North cell wall
    dudyS = (u[1:rows+1,1:cols+1]-u[0:cols,1:cols+1])/dy # du/dy on South cell wall
    
    dvdxE = (v[1:rows+1,2:]-v[1:rows+1,1:cols+1])/dx # dv/dx on East cell wall
    dvdxW = (v[1:rows+1,1:cols+1]-v[1:rows+1,0:cols])/dx # dv/dx on West cell wall
    dvdyN = (v[2:,1:cols+1]-v[1:rows+1,1:cols+1])/dy # dv/dy on North cell wall
    dvdyS = (v[1:rows+1,1:cols+1]-v[0:cols,1:cols+1])/dy # dv/dy on South cell wall

    # Complicated Differentials
    dvdxT = (v[2:,2:]-v[2:,0:cols])/(2*dx) # dvdx for j+1
    dvdxC = (v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx) # dvdx for j
    dvdxB = (v[0:rows,2:]-v[0:rows,0:cols])/(2*dx) # dvdx for j-1

    dudyR = (u[2:,2:]-u[0:rows,2:])/(2*dy) # dudy for i+1
    dudyC = (u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy) # dudy for i
    dudyL = (u[2:,0:cols]-u[0:rows,0:cols])/(2*dy) # dudy for i-1

    dvdxN = (dvdxT+dvdxC)/2 # dv/dx on North cell wall
    dvdxS = (dvdxC+dvdxB)/2 # dv/dx on South cell wall
    dudyE = (dudyR+dudyC)/2 # du/dy on East cell wall
    dudyW = (dudyC+dudyL)/2 # du/dy on West cell wall

    # Cell Wall Viscosities
    muW = (mu[1:rows+1,1:cols+1]+mu[1:rows+1,2:])/2 # Viscosity at West cell wall
    muE = (mu[1:rows+1,1:cols+1]+mu[1:rows+1,0:cols])/2
    muN = (mu[1:rows+1,1:cols+1]+mu[2:,1:cols+1])/2
    muS = (mu[1:rows+1,1:cols+1]+mu[0:cols,1:cols+1])/2

    # Convection Term
    conv_x = u[1:rows+1,1:cols+1]*dudx+v[1:rows+1,1:cols+1]*dudy
    conv_y = u[1:rows+1,1:cols+1]*dvdx+v[1:rows+1,1:cols+1]*dvdy
    # Diffusion Term
    diff_x = ( (muE*dudxE-muW*dudxW)*(2/dx) + (muN*(dudyN+dvdxN)-muS*(dudyS+dvdxS))*(1/dy) )/(Re*rho[1:rows+1,1:cols+1])
    diff_y = ( (muE*(dvdxE+dudyE)-muW*(dvdxW+dudyW))*(1/dx) + (muN*dvdyN-muS*dvdyS)*(2/dy) )/(Re*rho[1:rows+1,1:cols+1])
    # Gravitational Term
    grav_x = 0#gx/Fr**2
    grav_y = 0#gy/Fr**2
    # Surface Tension Term
    Bx = sigma*k[1:rows+1,1:cols+1]*delta[1:rows+1,1:cols+1]*dphidx*lc**2
    By = sigma*k[1:rows+1,1:cols+1]*delta[1:rows+1,1:cols+1]*dphidy*lc**2
    surf_x = Bx/(We*rho[1:rows+1,1:cols+1])
    surf_y = By/(We*rho[1:rows+1,1:cols+1])

    u_star = u.copy()
    v_star = v.copy()
    # u_star and v_star
    u_star[1:rows+1,1:cols+1] = u[1:rows+1,1:cols+1] + dt*(-conv_x+diff_x+grav_x+surf_x)
    v_star[1:rows+1,1:cols+1] = v[1:rows+1,1:cols+1] + dt*(-conv_y+diff_y+grav_y+surf_y)

    print(f"dt*conv_x: {np.amax(dt*conv_x)}")
    print(f"dt*diff_x: {np.amax(dt*diff_x)}")
    print(f"dt*grav_x: {np.amax(dt*grav_x)}")
    print(f"dt*surf_x: {np.amax(dt*surf_x)}")
    print(f"dt*conv_y: {np.amax(dt*conv_y)}")
    print(f"dt*diff_y: {np.amax(dt*diff_y)}")
    print(f"dt*grav_y: {np.amax(dt*grav_y)}")
    print(f"dt*surf_y: {np.amax(dt*surf_y)}")

    space.u_star = u_star.copy()
    space.v_star = v_star.copy()