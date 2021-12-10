"""
Module to Handle Magnetic Fields and some other general imports
"""
import numpy as np 
import matplotlib.pyplot as plt 

class Quadrapole():
    """
    Implements Quadrapole magnetic fields  
    Attributes
    ----------
    G0 - Peak field gradient 
    L - Geometric length of QP 
    R0 - Bore radius of QP
    cis - Enge's coef
    zi - z-coordinate of Start of QP 
    zf - z-coordinate of End of QP
    zc - z-coordinate of Center of QP

    Methods
    -------
    field(x,y,z) 
        Returns Bx, By, Bz  components of QP Magnetic Field
    """
    def __init__(self, G0=20, L=1, R0=1, zi=0) : 
        
        self.G0 = G0 
        self.L = L 
        self.R0 = R0
        self.zi = zi
        self.zc = zi + L/2
        self.zf = zi + L 

       
        self.lamda = 2*self.R0

    def field(self, x, y, z): 
        c0 = 0.1122
        c1 = 6.2671 
        c2 = -1.4982
        c3 = 3.5882
        c4 = -2.1209 
        c5 = 1.723
        L = self.L/2 
        lamda = self.lamda
        Nr = self.G0
        z0 = self.zc
        z=z-z0

        if z<0:
            k = (z+L)/lamda
        else:
            k= (z-L)/lamda

        Dr = c0 + c1*k + c2*k*k + c3*k**3 + c4*k**4 + c5*k**5
        Dr = 1 + np.exp(Dr)
        G = Nr/Dr
        dDr = (Dr-1)*(c1 + 2*c2*k + 3*c3*k*k + 4*c4*k**3 + 5*c5*k**4)/lamda

        Bx = G*y 
        By = G*x 
        Bz = -G*x*y*dDr/Dr

        B = np.array([Bx, By, Bz])
        return B

    def __str__(self): 
        return "Quadrupole"

class Drift():
    """
    Implements Drift Space  
    Attributes
    ---------- 
    L - Geometric length of drift space 
    zi - z-coordinate of Start of Drift space 
    zf - z-coordinate of End of Drift space
    zc - z-coordinate of Center of Drift space

    Methods
    -------
    field(x,y,z) 
        Returns Bx, By, Bz  components of Drift space i.e 0,0,0
    """
    def __init__(self, L=1, zi=0, R0=1) :
        self.L = L 
        self.zc = zi + L/2
        self.zi = zi 
        self.zf = zi + L
        self.R0 = R0

    def field(self, x, y, z):
        Bx = 0
        By = 0
        Bz = 0
        B = np.array([Bx, By, Bz])
        return B
    
    def __str__(self): 
        return "Drift"

class Mirror():
    """
    Implements Quadrapole magnetic fields  
    Attributes
    ----------
    B0 - Peak field
    L -  Half Geometric length of Magnetic Mirror
    R0 - Bore radius of QP
    zi - z-coordinate of Start of Mirror
    zf - z-coordinate of End of Mirror
    Zc - z-coordinate of Center of Mirror

    Methods
    -------
    field(x,y,z) 
        Returns Bx, By, Bz  components of Magnetic Mirror Field
    """
    def __init__(self, B0, L, zi):
        self.B0 = B0
        self.L = L
        self.zi = zi 
        self.zc = zi + L/2
        self.zf = zi + L
        
    def field(self,x,y,z):
        B0 = self. B0 
        L = self.L/2 
        a = 2*L
        Bx = np.array(x*B0*z, dtype=np.float64)
        Bx = -1*Bx/L**2
        By = np.array(y*B0*z, dtype=np.float64)
        By = -1*By/L**2
        Bz = np.array(B0*(1+z**2/L**2), dtype=np.float64)
        B = (Bx,By,Bz)
        return B
       
    def __str__(self): 
        return "Mirror"

class Trajactory():
    """
    Contains the trajectory of the particle

    Attributes
    ----------
    t - Array containing time slices 
    X - Array containing position velocity vector for each time slice  
    

    Methods
    -------
    plotXt(sol)
        Plots variation of x,y,z with time t     
    
    plotXY(sol) 
        Plots variation of x,y,z with time t 
      
    """
    def __init__(self, t, X) :
        self.t = t 
        self.X = X 

    def plotXt(self,sol= None):
        if sol is None: 
            t = self.t 
            x, y, z = self.x, self.y, self.z
        else:
            t,X = sol
            x = X[:,0]
            y = X[:,1]
            z = X[:,2]
        
        fig, (ax1,ax2,ax3) = plt.subplots(3,1)
        fig.set_figwidth(10)
        ax1.set_title("Position vs Time")
        ax1.plot(t,x)
        ax1.set_ylabel("x")
        ax2.plot(t,y)
        ax2.set_ylabel("y")
        ax3.plot(t,z)
        ax3.set_ylabel("z")
        ax2.set_xlabel("t")
  
    def plotXY(self,sol=None):
        if sol is None:
            x, y, z = self.x, self.y, self.z
        else:
            t,X = sol
            x = X[:,0]
            y = X[:,1]
            z = X[:,2]

        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.set_figwidth(10)
        ax1.set_title("Motion in xy plane")
        ax1.plot(x,y)
        ax1.set_ylabel("y")
        ax1.set_xlabel("x")
        ax2.plot(x,z)
        ax2.set_title("Motion in xz plane")
        ax2.set_ylabel("z")
        ax2.set_xlabel("x")


# -------------------------------- To be Completed ------------------------------------------- #
# TODO: 
class Solenoid():
    """
    Implements Quadrapole magnetic fields  
    Attributes
    ----------
    B0 - Peak field
    L -  Half Geometric length of Solenoid
    zi - z-coordinate of Start of Mirror
    zf - z-coordinate of End of Mirror
    Zc - z-coordinate of Center of Mirror

    Methods
    -------
    field(x,y,z) 
        Returns Bx, By, Bz  components of Magnetic Mirror Field
    """
    def __init__(self, B0, L, zi): 
        self.B0 = B0
        self.L = L
        self.zi = zi 
        self.zc = zi + L/2
        self.zf = zi + L 
        
    def field(self,x,y,z):
        pass
    
    def __str__(self): 
        return "Solenoid"

class FieldPlotter():
    """
    Plots the Magnetic fields
    """
    def __init__(self, MagneticElement, N=100) -> None:
        self.L = MagneticElement.L   
        self.a = 2*self.L
        self.field = MagneticElement.field 
        self.N = N 

    def createXYMesh(self,z0=1):
        a = self.a 
        x = np.linspace(-a,a,self.N)
        y = np.linspace(-a,a,self.N)
        x_mesh,y_mesh = np.meshgrid(x,y)
        z = np.ones(len(x_mesh))
        z = z-1+z0
        xy_mesh = (x_mesh,y_mesh,z)
        return xy_mesh

    def createXZMesh(self,y0=1):
        a = self.a 
        x = np.linspace(-a,a,self.N)
        z = np.linspace(-a,a,self.N)
        x_mesh,z_mesh = np.meshgrid(x,z)
        y = np.ones(len(x_mesh))
        y = y-1+y0
        xz_mesh = (x_mesh,y,z_mesh)
        return xz_mesh
        
    def createYZMesh(self,x0=1):
        a = self.a 
        y = np.linspace(-a,a,self.N)
        z = np.linspace(-a,a,self.N)
        y_mesh,z_mesh = np.meshgrid(y,z)
        x = np.ones(len(y_mesh))
        x = x-1+x0
        yz_mesh = (x,y_mesh,z_mesh)
        return yz_mesh

    def plotMageticFieldXY(self, z0=1, magneticField = None):
        if magneticField is None: magneticField = self.field
        x,y,z = self.createXYMesh(z0)
        B = magneticField(x,y,z)
        plt.title("Magnetic Field in x-y plane at z = "+ str(z0))
        magnitude = np.sqrt(B[0]**2 + B[1]**2)
        lw = 3*magnitude/np.amax(magnitude)
        plt.streamplot(x,y,B[0],B[1],linewidth=lw)

    
    def plotMageticFieldXZ(self,y0=1, magneticField = None):
        if magneticField is None: magneticField = self.field
        x,y,z = self.createXZMesh(y0)
        B = magneticField(x,y,z)
        plt.title("Magnetic Field in x-z plane at y = "+ str(y0))
        magnitude = np.sqrt(B[0]**2 + B[2]**2)
        lw = 3*magnitude/np.amax(magnitude)
        plt.streamplot(x,z,B[0],B[2],linewidth=lw)
    
    def plotMageticFieldYZ(self,x0=1, magneticField = None):
        if magneticField is None: magneticField = self.field
        x,y,z = self.createYZMesh(x0)
        B = magneticField(x,y,z)
        plt.title("Magnetic Field in y-z plane at x = "+ str(x0))
        magnitude = np.sqrt(B[1]**2 + B[2]**2)
        lw = 3*magnitude/np.amax(magnitude)
        plt.streamplot(y,z,B[1],B[2],linewidth=lw)

class OpticalLayout():
    def __init__(self) -> None:
        self.MAGNETICELTSNAMEDICT = {"Drift":3,"Quadrapole":5} 
        self.MAGNETICELTSCLASSDICT = {3: Drift, 5:Quadrapole}

    def readLayout(self, fname="Optics_layout.txt"):
        with open(fname) as f:
            pass 


